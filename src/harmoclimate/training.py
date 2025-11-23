"""Linear harmonic model training utilities for HarmoClimate."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
from scipy.stats import norm, skewnorm

from .config import AUTHOR_NAME, COUNTRY_CODE, MODEL_VERSION
from .core import SOLAR_YEAR_DAYS, prepare_dataset


# ----------------------------- Global settings -----------------------------

FINAL_TRAINING_PERIOD_LABEL = "1999–2025"
RIDGE_LAMBDA_DEFAULT = 0.0
DEFAULT_DAYS_INCLUSIVE_MAX = int(math.floor(SOLAR_YEAR_DAYS))

# Target variable constants
TARGET_TEMP = "T"
TARGET_HUMIDITY = "Q"
TARGET_PRESSURE = "P"

UNIT_TEMP = "degC"
UNIT_HUMIDITY = "kg/kg"
UNIT_PRESSURE = "hPa"

# ----------------------------- Data classes -------------------------------


@dataclass
class ErrorMetrics:
    """Collection of scalar error metrics for a predicted variable."""

    mae: float
    bias: float
    err_p05: float
    err_p95: float

    def as_dict(self) -> dict[str, float]:
        return {
            "mae": self.mae,
            "bias": self.bias,
            "err_p05": self.err_p05,
            "err_p95": self.err_p95,
        }


@dataclass
class YearlyDesignStats:
    """Per-year sufficient statistics for a given target variable."""

    year: int
    X: np.ndarray
    y: np.ndarray
    S: np.ndarray
    b: np.ndarray
    n: int
    utc_day_index: np.ndarray
    utc_hour: np.ndarray
    params_meta: list[dict[str, int]]


@dataclass
class YearlyValidationMetrics:
    """Per-year validation diagnostics for the leave-one-year-out protocol."""

    year: int
    mse_model: float
    mse_ref: float
    rmse: float
    skill: float
    n: int

    def as_dict(self) -> dict[str, object]:
        return {
            "year": int(self.year),
            "mse_model": float(self.mse_model),
            "mse_ref": float(self.mse_ref),
            "rmse": float(self.rmse),
            "skill": float(self.skill),
            "n": int(self.n),
        }


@dataclass
class LeaveOneYearOutReport:
    """Summary of leave-one-year-out validation for a given target."""

    years: list[YearlyValidationMetrics]
    global_rmse: float
    global_skill: float
    total_observations: int
    hyperparameters: dict[str, object]
    ridge_lambda: float
    final_training_period: str
    protocol: str = "leave_one_year_out"
    evaluation_time_base: str = "UTC"
    model_time_base: str = "solar"
    baseline: str = "climatology_mean per (utc_day, utc_hour), LOYO"

    def to_payload(self) -> dict[str, object]:
        return {
            "protocol": self.protocol,
            "final_training_period": self.final_training_period,
            "hyperparameters": dict(self.hyperparameters),
            "global": {
                "rmse": float(self.global_rmse),
                "skill": float(self.global_skill),
            },
            "years_validated": [entry.as_dict() for entry in self.years],
        }


@dataclass
class ParameterLayout:
    """Describes how a slice of the coefficient vector maps to a model parameter."""

    name: str
    role: str
    n_annual: int
    start: int
    length: int

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "role": self.role,
            "n_annual": self.n_annual,
            "start": self.start,
            "length": self.length,
        }


@dataclass
class LinearModelFit:
    """Fitted linear harmonic model for a single target variable."""

    target_variable: str
    target_unit: str
    coefficients: np.ndarray
    params_layout: list[ParameterLayout]
    n_diurnal: int
    default_n_annual: int
    annual_per_param: dict[str, int]
    metrics: ErrorMetrics
    validation: Optional[LeaveOneYearOutReport] = None

    def coefficients_list(self) -> list[float]:
        return [float(v) for v in self.coefficients]


@dataclass
class TrainingResult:
    """Bundle capturing linear models for temperature, specific humidity, and pressure."""

    temperature_model: LinearModelFit
    specific_humidity_model: LinearModelFit
    pressure_model: LinearModelFit


@dataclass
class StochasticModelFit:
    """Fitted VAR(p) stochastic model for residual dependencies between P, Q, T."""

    coefficient_matrix: np.ndarray  # (3*len(lags)) x 3 coefficient matrix (A)
    covariance_matrix: np.ndarray  # 3x3 noise covariance matrix (Sigma)
    variables: list[str]  # Variable names in order: ["P", "Q", "T"]
    skew_normal_params: Dict[str, Dict[str, List[float]]]  # Skew normal harmonic coeffs for each variable
    model_type: str = "VAR"
    lags: list[int] = (1, 2)

    def coefficient_matrix_list(self) -> list[list[float]]:
        """Convert coefficient matrix to nested list for JSON serialization."""
        return [[float(v) for v in row] for row in self.coefficient_matrix]

    def covariance_matrix_list(self) -> list[list[float]]:
        """Convert covariance matrix to nested list for JSON serialization."""
        return [[float(v) for v in row] for row in self.covariance_matrix]


# ----------------------------- Utilities ----------------------------------


def _compute_error_metrics(errors: np.ndarray) -> ErrorMetrics:
    """Compute basic error diagnostics for residuals."""

    valid = errors[np.isfinite(errors)]
    if valid.size == 0:
        return ErrorMetrics(mae=math.nan, bias=math.nan, err_p05=math.nan, err_p95=math.nan)
    mae = float(np.mean(np.abs(valid)))
    bias = float(np.mean(valid))
    err_p05 = float(np.quantile(valid, 0.05))
    err_p95 = float(np.quantile(valid, 0.95))
    return ErrorMetrics(mae=mae, bias=bias, err_p05=err_p05, err_p95=err_p95)


def solve_normal_equations(S: np.ndarray, b: np.ndarray, ridge_lambda: float) -> np.ndarray:
    """Solve (S + λI) β = b with a numerically stable fallback."""

    if S.ndim != 2 or S.shape[0] != S.shape[1]:
        raise ValueError("Matrix S must be square for normal equation solving.")
    if b.ndim != 1 or b.shape[0] != S.shape[0]:
        raise ValueError("Vector b must align with S dimensions.")

    S_reg = S.astype(float, copy=True)
    if ridge_lambda > 0.0:
        S_reg = S_reg + ridge_lambda * np.eye(S_reg.shape[0], dtype=float)

    try:
        return np.linalg.solve(S_reg, b)
    except np.linalg.LinAlgError:
        solution, *_ = np.linalg.lstsq(S_reg, b, rcond=None)
        return solution


def prepare_training_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Attach solar descriptors, derived fields, and calendar metadata."""

    working = prepare_dataset(df)
    working["year"] = working["DT_UTC"].dt.year.astype("Int64")
    working = working.dropna(subset=["year"]).copy()
    working["year"] = working["year"].astype(int)
    return working


# ----------------------------- Design matrix ------------------------------


def build_annual_basis(day: np.ndarray, n_annual: int) -> np.ndarray:
    """Construct the annual harmonic basis (constant + cos/sin pairs)."""

    omega = 2.0 * math.pi / SOLAR_YEAR_DAYS
    cols = [np.ones_like(day)]
    for k in range(1, n_annual + 1):
        cols.append(np.cos(k * omega * day))
        cols.append(np.sin(k * omega * day))
    return np.column_stack(cols)


def _role_for_parameter(name: str) -> str:
    """Determine the role of a parameter based on its name."""
    if name == "c0":
        return "offset"
    if name.startswith("a"):
        return f"diurnal_cos_{int(name[1:])}"
    if name.startswith("b"):
        return f"diurnal_sin_{int(name[1:])}"
    raise ValueError(f"Unknown parameter name '{name}'")


def _prepare_time_indices(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Compute UTC day index and hour from the dataframe's DT_UTC column.

    Handles leap years and ensures timestamps are hourly.
    """
    utc_requirement_msg = (
        "LOYO evaluation runs in UTC and requires hourly UTC timestamps (`DT_UTC` at whole hours)."
    )
    if "DT_UTC" not in df.columns:
        raise KeyError(utc_requirement_msg)
    
    utc_series = pd.to_datetime(df["DT_UTC"], utc=True, errors="coerce")
    if utc_series.isna().any() or not utc_series.equals(utc_series.dt.floor("h")):
        raise ValueError(utc_requirement_msg)
    
    month = np.asarray(utc_series.dt.month, dtype=int)
    day_of_month = np.asarray(utc_series.dt.day, dtype=int)
    hour_utc = np.asarray(utc_series.dt.hour, dtype=int)
    day_of_year = np.asarray(utc_series.dt.dayofyear, dtype=int)
    is_leap_year = np.asarray(utc_series.dt.is_leap_year, dtype=bool)

    utc_day_index = day_of_year.astype(int, copy=True)
    leap_adjust_mask = is_leap_year & (utc_day_index > 59)
    utc_day_index[leap_adjust_mask] -= 1
    feb29_mask = (month == 2) & (day_of_month == 29)
    utc_day_index[feb29_mask] = -1
    utc_hour = hour_utc.astype(int, copy=True)
    
    return utc_day_index, utc_hour


def build_global_linear_matrix(
    df: pd.DataFrame,
    *,
    n_diurnal: int,
    annual_per_param: Dict[str, int],
    default_n_annual: int,
    target: str,
) -> tuple[
    np.ndarray,
    np.ndarray,
    list[dict[str, int]],
    np.ndarray,
    np.ndarray,
]:
    """Construct the design matrix for the factorized linear model.

    UTC no-leap calendar for evaluation keys: Feb 29 samples are excluded from
    `(utc_day_index, utc_hour)`.

    Returns:
        (
            X,
            y,
            params_meta,
            utc_day_index,
            utc_hour,
        )
    """

    required_cols = {"yday_frac_solar", "hour_solar", target}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns for training: {sorted(missing)}")

    local = df.copy()
    local = local.dropna(subset=["yday_frac_solar", "hour_solar", target])

    local["day"] = local["yday_frac_solar"].astype(float)
    local["hour"] = np.mod(local["hour_solar"].astype(float), 24.0)

    day = local["day"].to_numpy(dtype=float)
    hour = local["hour"].to_numpy(dtype=float)
    y = local[target].to_numpy(dtype=float)

    utc_day_index, utc_hour = _prepare_time_indices(local)

    omega = 2.0 * math.pi / 24.0
    param_names: list[str] = ["c0"]
    for m in range(1, n_diurnal + 1):
        param_names.extend([f"a{m}", f"b{m}"])

    X_blocks: list[np.ndarray] = []
    params_meta: list[dict[str, int]] = []
    start_idx = 0

    for name in param_names:
        n_annual_param = annual_per_param.get(name, default_n_annual)
        annual_basis = build_annual_basis(day, n_annual_param)

        if name == "c0":
            block = annual_basis
        elif name.startswith("a"):
            harmonic = int(name[1:])
            block = annual_basis * np.cos(harmonic * omega * hour)[:, None]
        elif name.startswith("b"):
            harmonic = int(name[1:])
            block = annual_basis * np.sin(harmonic * omega * hour)[:, None]
        else:
            raise ValueError(f"Unhandled parameter name '{name}'")

        X_blocks.append(block)
        length = annual_basis.shape[1]
        params_meta.append(
            {"name": name, "n_annual": n_annual_param, "start": start_idx, "length": length}
        )
        start_idx += length

    X = np.concatenate(X_blocks, axis=1)
    return X, y, params_meta, utc_day_index, utc_hour


def _build_layout(meta: Iterable[dict[str, int]]) -> list[ParameterLayout]:
    """Convert raw parameter metadata into a structured layout list."""
    layout: list[ParameterLayout] = []
    for entry in meta:
        name = entry["name"]
        role = _role_for_parameter(name)
        layout.append(
            ParameterLayout(
                name=name,
                role=role,
                n_annual=int(entry["n_annual"]),
                start=int(entry["start"]),
                length=int(entry["length"]),
            )
        )
    return layout


# ----------------------------- Training helpers ---------------------------


def compute_sufficient_stats(
    df: pd.DataFrame,
    *,
    target: str,
    n_diurnal: int,
    default_n_annual: int,
    annual_per_param: dict[str, int] | None,
) -> list[YearlyDesignStats]:
    """Compute per-year design matrices and sufficient statistics."""

    if target not in df.columns:
        raise KeyError(f"Target column '{target}' is missing from the training frame.")
    if "year" not in df.columns:
        raise KeyError("Training frame is missing the 'year' column; call prepare_training_frame.")

    annual_per_param = dict(annual_per_param or {})
    stats: list[YearlyDesignStats] = []

    for year, group in df.groupby("year"):
        if group.empty:
            continue

        (
            X,
            y,
            meta,
            utc_day_index,
            utc_hour,
        ) = build_global_linear_matrix(
            group,
            n_diurnal=n_diurnal,
            annual_per_param=annual_per_param,
            default_n_annual=default_n_annual,
            target=target,
        )
        if y.size == 0 or X.size == 0:
            continue

        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)

        S = X.T @ X
        b = X.T @ y

        stats.append(
            YearlyDesignStats(
                year=int(year),
                X=X,
                y=y,
                S=S,
                b=b,
                n=int(y.size),
                utc_day_index=utc_day_index.astype(int, copy=False),
                utc_hour=utc_hour.astype(int, copy=False),
                params_meta=[dict(entry) for entry in meta],
            )
        )

    return stats


def fit_from_stats(
    stats: List[YearlyDesignStats],
    ridge_lambda: float,
) -> tuple[np.ndarray, ErrorMetrics, list[ParameterLayout]]:
    """Fit coefficients from sufficient statistics and return diagnostics."""

    if not stats:
        raise ValueError("No sufficient statistics supplied for fitting.")

    feature_dim = stats[0].S.shape[0]
    S_total = np.zeros((feature_dim, feature_dim), dtype=float)
    b_total = np.zeros(feature_dim, dtype=float)

    for entry in stats:
        if entry.S.shape != (feature_dim, feature_dim):
            raise ValueError("Inconsistent feature dimension across yearly statistics.")
        S_total += entry.S
        b_total += entry.b

    coefficients = solve_normal_equations(S_total, b_total, ridge_lambda)

    residuals: list[np.ndarray] = []
    for entry in stats:
        y_pred = entry.X @ coefficients
        valid_mask = np.isfinite(entry.y) & np.isfinite(y_pred)
        if np.any(valid_mask):
            residuals.append(entry.y[valid_mask] - y_pred[valid_mask])

    all_residuals = np.concatenate(residuals) if residuals else np.empty(0, dtype=float)
    metrics = _compute_error_metrics(all_residuals)
    layout = _build_layout(stats[0].params_meta)
    return coefficients.astype(float), metrics, layout


def _is_prepared(df: pd.DataFrame) -> bool:
    """Check if the dataframe has all required columns for training."""
    required = {
        "yday_frac_solar",
        "hour_solar",
        "year",
        "Q",
    }
    return required.issubset(df.columns)


def _train_target(
    df: pd.DataFrame,
    *,
    target_variable: str,
    target_unit: str,
    n_diurnal: int,
    default_n_annual: int,
    annual_per_param: dict[str, int],
    ridge_lambda: float,
) -> LinearModelFit:
    """Train a single target variable using the provided configuration."""
    stats = compute_sufficient_stats(
        df,
        target=target_variable,
        n_diurnal=n_diurnal,
        default_n_annual=default_n_annual,
        annual_per_param=annual_per_param,
    )
    if not stats:
        raise ValueError(f"No samples available to train target '{target_variable}'.")

    coefficients, metrics, layout = fit_from_stats(stats, ridge_lambda)

    return LinearModelFit(
        target_variable=target_variable,
        target_unit=target_unit,
        coefficients=coefficients,
        params_layout=layout,
        n_diurnal=n_diurnal,
        default_n_annual=default_n_annual,
        annual_per_param=dict(annual_per_param),
        metrics=metrics,
        validation=None,
    )


def train_models(
    df: pd.DataFrame,
    *,
    n_diurnal: int = 3,
    default_n_annual: int = 3,
    annual_per_param: dict[str, int] | None = None,
    ridge_lambda: float = RIDGE_LAMBDA_DEFAULT,
) -> TrainingResult:
    """Model is trained on solar features; no external evaluation is performed here."""

    annual_per_param = dict(annual_per_param or {})

    working = df if _is_prepared(df) else prepare_training_frame(df)

    temperature_model = _train_target(
        working,
        target_variable=TARGET_TEMP,
        target_unit=UNIT_TEMP,
        n_diurnal=n_diurnal,
        default_n_annual=default_n_annual,
        annual_per_param=annual_per_param,
        ridge_lambda=ridge_lambda,
    )
    specific_humidity_model = _train_target(
        working,
        target_variable=TARGET_HUMIDITY,
        target_unit=UNIT_HUMIDITY,
        n_diurnal=n_diurnal,
        default_n_annual=default_n_annual,
        annual_per_param=annual_per_param,
        ridge_lambda=ridge_lambda,
    )
    pressure_model = _train_target(
        working,
        target_variable=TARGET_PRESSURE,
        target_unit=UNIT_PRESSURE,
        n_diurnal=n_diurnal,
        default_n_annual=default_n_annual,
        annual_per_param=annual_per_param,
        ridge_lambda=ridge_lambda,
    )

    return TrainingResult(
        temperature_model=temperature_model,
        specific_humidity_model=specific_humidity_model,
        pressure_model=pressure_model,
    )


def _gamma1_from_delta(delta: np.ndarray) -> np.ndarray:
    """
    Theoretical skewness of a skew-normal as a function of delta.
    delta: array-like, values in (-1, 1).
    """
    delta = np.asarray(delta, dtype=float)
    delta2 = delta * delta

    m1 = delta * math.sqrt(2.0 / math.pi)
    v1 = 1.0 - 2.0 * delta2 / math.pi

    num = (4.0 - math.pi) / 2.0 * (m1**3)
    denom = np.power(v1, 1.5)
    return num / denom


def invert_skewness_to_delta(g1: np.ndarray,
                             max_iter: int = 30,
                             tol: float = 1e-8) -> np.ndarray:
    """
    Invert empirical skewness g1 to delta for a skew-normal distribution.
    Vectorized over all entries of g1.
    """
    g1 = np.asarray(g1, dtype=float)

    # Initial guess: small skewness -> small delta
    delta = np.sign(g1) * np.power(np.abs(g1) / 2.0 + 1e-12, 1.0 / 3.0)
    delta = np.clip(delta, -0.95, 0.95)

    for _ in range(max_iter):
        g_est = _gamma1_from_delta(delta)
        diff = g_est - g1

        if np.all(np.abs(diff) < tol):
            break

        # Numerical derivative dγ/dδ via finite differences (vectorized)
        h = 1e-5
        d_plus = np.clip(delta + h, -0.99, 0.99)
        d_minus = np.clip(delta - h, -0.99, 0.99)

        g_plus = _gamma1_from_delta(d_plus)
        g_minus = _gamma1_from_delta(d_minus)
        deriv = (g_plus - g_minus) / (2.0 * h)

        deriv = np.where(np.abs(deriv) < 1e-10, 1e-10, deriv)

        delta_new = delta - diff / deriv
        delta = np.clip(delta_new, -0.99, 0.99)

    return delta


def delta_to_a(delta: np.ndarray) -> np.ndarray:
    """
    Convert delta to shape parameter a.
    """
    delta = np.asarray(delta, dtype=float)
    return delta / np.sqrt(1.0 - delta * delta)


def train_stochastic_model(
    training_result: TrainingResult,
    df: pd.DataFrame,
    lags: Sequence[int] = (1, 6, 12),
) -> StochasticModelFit:
    """Train a VAR(p) stochastic model on residuals from the harmonic models.
    
    This function:
    1. Predicts baseline values using the harmonic models for P, Q, T
    2. Computes residuals: observed - predicted
    3. Fits annual harmonic envelopes for skew-normal parameters (shape a, sigma)
       to the residuals, enforcing zero mean.
    4. Transforms residuals to standard normal Z-scores.
    5. Trains a sparse VAR model on Z-scores.
    
    Args:
        training_result: Fitted harmonic models for P, Q, T
        df: Prepared dataframe with solar time features and observed values
        lags: List of lags to include in the model (e.g. [1, 2, 24])
        
    Returns:
        StochasticModelFit with VAR coefficient matrix A and covariance Sigma
    """
    from .core import compute_solar_time, SOLAR_YEAR_DAYS, skew_params_from_a_sigma
    
    # Ensure dataframe is prepared
    working = df if _is_prepared(df) else prepare_training_frame(df)
    
    # Get longitude for solar time computation
    if "LON" not in working.columns:
        raise KeyError("Column 'LON' required for solar time computation")
    lon = working["LON"].iloc[0]
    
    # Predict harmonic baselines for each variable
    def predict_harmonic(model: LinearModelFit, timestamps: pd.Series) -> np.ndarray:
        """Predict values using the trained harmonic model."""
        # Compute solar time features
        lon_series = pd.Series(data=lon, index=timestamps.index)
        solar = compute_solar_time(timestamps, lon_series)
        day = solar["yday_frac_solar"].to_numpy(dtype=float)
        hour = solar["hour_solar"].to_numpy(dtype=float)
        
        omega_day = 2.0 * math.pi / SOLAR_YEAR_DAYS
        omega_hour = 2.0 * math.pi / 24.0
        
        # Build design matrix from model layout
        y_pred = np.zeros(len(timestamps))
        coeffs = model.coefficients
        
        for param in model.params_layout:
            # Build annual basis
            annual_basis_cols = [np.ones_like(day)]
            for k in range(1, param.n_annual + 1):
                annual_basis_cols.append(np.cos(k * omega_day * day))
                annual_basis_cols.append(np.sin(k * omega_day * day))
            annual_basis = np.column_stack(annual_basis_cols)
            
            # Apply diurnal modulation
            if param.name == "c0":
                block = annual_basis
            elif param.name.startswith("a"):
                h = int(param.name[1:])
                block = annual_basis * np.cos(h * omega_hour * hour)[:, None]
            elif param.name.startswith("b"):
                h = int(param.name[1:])
                block = annual_basis * np.sin(h * omega_hour * hour)[:, None]
            else:
                raise ValueError(f"Unknown parameter name: {param.name}")
            
            # Apply coefficients
            block_coeffs = coeffs[param.start : param.start + param.length]
            y_pred += block @ block_coeffs
        
        return y_pred
    
    # Compute residuals for each variable
    residuals_df = pd.DataFrame(index=working.index)
    
    # Pressure
    baseline_P = predict_harmonic(training_result.pressure_model, working["DT_UTC"])
    residuals_df["P"] = working["P"].to_numpy() - baseline_P
    
    # Specific humidity
    baseline_Q = predict_harmonic(training_result.specific_humidity_model, working["DT_UTC"])
    residuals_df["Q"] = working["Q"].to_numpy() - baseline_Q
    
    # Temperature
    baseline_T = predict_harmonic(training_result.temperature_model, working["DT_UTC"])
    residuals_df["T"] = working["T"].to_numpy() - baseline_T
    
    # Drop NaNs
    residuals_df = residuals_df.dropna()
    
    # Extract residual matrix in order [P, Q, T]
    variable_order = ["P", "Q", "T"]
    skew_params_coeffs = {}
    z_scores_list = []
    
    # Solar day for envelope fitting
    # Align day_solar with residuals (after dropna)
    day_solar = working.loc[residuals_df.index, "yday_frac_solar"].to_numpy(dtype=float)
    
    day_index = np.floor(day_solar).astype(int)
    day_index = np.clip(day_index, 1, 365)
    
    n_days = 365
    MIN_POINTS = 50
    
    for var in variable_order:
        residuals = residuals_df[var].to_numpy(dtype=float)
        
        # 1. Daily empirical estimates
        sigma_raw = np.full(n_days, np.nan, dtype=float)
        g1_raw = np.full(n_days, np.nan, dtype=float)
        
        for d in range(1, n_days + 1):
            mask = (day_index == d)
            subset = residuals[mask]
            subset = subset[np.isfinite(subset)]
            
            if subset.size >= MIN_POINTS:
                # Enforce local zero-mean
                m = np.mean(subset)
                subset_centered = subset - m
                
                m2 = np.mean(subset_centered**2)
                m3 = np.mean(subset_centered**3)
                
                eps = 1e-12
                sigma_d = math.sqrt(m2 + eps)
                g1_d = m3 / (sigma_d**3 + eps)
                
                idx = d - 1
                sigma_raw[idx] = sigma_d
                g1_raw[idx] = g1_d
        
        # 2. Convert skewness to shape a
        valid = np.isfinite(g1_raw) & np.isfinite(sigma_raw)
        g1_daily = g1_raw.copy()
        g1_daily[~valid] = 0.0 # Neutral value
        
        delta_daily = invert_skewness_to_delta(g1_daily)
        a_raw = delta_to_a(delta_daily)
        
        # 3. Fit annual harmonic envelopes
        # We use 2 harmonics for envelopes to capture seasonality
        n_envelope_harmonics = 2
        
        grid_days = np.arange(1, 366, dtype=float)
        basis = build_annual_basis(grid_days, n_annual=n_envelope_harmonics)
        
        # Restrict to valid days
        mask_valid = np.isfinite(a_raw) & np.isfinite(sigma_raw) & valid
        if np.sum(mask_valid) < 10:
             # Fallback if too few valid days
             coeffs_a = np.zeros(1 + 2 * n_envelope_harmonics)
             coeffs_sigma = np.zeros(1 + 2 * n_envelope_harmonics)
             coeffs_sigma[0] = np.nanstd(residuals) # Constant term
        else:
            basis_valid = basis[mask_valid]
            a_valid = a_raw[mask_valid]
            sigma_valid = sigma_raw[mask_valid]
            
            coeffs_a, _, _, _ = np.linalg.lstsq(basis_valid, a_valid, rcond=None)
            coeffs_sigma, _, _, _ = np.linalg.lstsq(basis_valid, sigma_valid, rcond=None)
            
        skew_params_coeffs[var] = {
            "a": [float(c) for c in coeffs_a],
            "sigma": [float(c) for c in coeffs_sigma],
        }
        
        # 4. Transform to Z-scores
        basis_full = build_annual_basis(day_solar, n_annual=n_envelope_harmonics)
        a_ts = basis_full @ coeffs_a
        sigma_ts = basis_full @ coeffs_sigma
        
        loc_ts, scale_ts = skew_params_from_a_sigma(a_ts, sigma_ts)
        
        u = skewnorm.cdf(residuals, a_ts, loc=loc_ts, scale=scale_ts)
        eps = 1e-6
        u = np.clip(u, eps, 1.0 - eps)
        z = norm.ppf(u)
        z_scores_list.append(z)

    # Stack Z-scores into matrix (n_samples, 3)
    X = np.column_stack(z_scores_list)
    
    # Prepare lag matrices for sparse VAR
    n_samples, n_vars = X.shape
    max_lag = max(lags) if lags else 0
    
    if n_samples <= max_lag:
        raise ValueError(f"Not enough samples ({n_samples}) to train VAR with max lag {max_lag}.")

    X_target = X[max_lag:]

    lagged_blocks = []
    for lag in lags:
        # Shift X by lag.
        block = X[max_lag - lag : n_samples - lag]
        lagged_blocks.append(block)

    X_lag = np.hstack(lagged_blocks)
    
    # Solve for A using least squares: A = lstsq(X_lag, X_target)
    A, _, _, _ = np.linalg.lstsq(X_lag, X_target, rcond=None)
    
    # Compute residuals and noise covariance
    X_pred = X_lag @ A
    E = X_target - X_pred
    Sigma = np.cov(E.T)
    
    return StochasticModelFit(
        coefficient_matrix=A.astype(float),
        covariance_matrix=Sigma.astype(float),
        variables=variable_order,
        skew_normal_params=skew_params_coeffs,
        model_type="VAR_sparse_harmonic_skew",
        lags=list(lags),
    )


# ----------------------------- Payload export -----------------------------


def build_parameter_payload(
    model: LinearModelFit,
    metadata: dict[str, object],
    generation_date_utc: str,
) -> dict[str, object]:
    """Build the JSON-friendly payload for a single fitted model."""

    payload = {
        "metadata": {
            "version": MODEL_VERSION,
            "generated_at_utc": generation_date_utc,
            "country_code": COUNTRY_CODE,
            "author": AUTHOR_NAME,
            "target_variable": model.target_variable,
            "target_unit": model.target_unit,
            **metadata,
        },
        "model": {
            "n_diurnal": model.n_diurnal,
            "params_layout": [entry.as_dict() for entry in model.params_layout],
            "coefficients": model.coefficients_list(),
        },
    }
    payload["metadata"]["error_envelope"] = {
        "mae": model.metrics.mae,
        "bias": model.metrics.bias,
        "p05": model.metrics.err_p05,
        "p95": model.metrics.err_p95,
    }
    payload["metadata"]["time_basis"] = {
        "type": "solar",
        "days": SOLAR_YEAR_DAYS,
        "calendar": "no-leap",
    }
    if model.validation is not None:
        payload["metadata"]["training_loyo_rmse"] = model.validation.global_rmse
        payload["metadata"]["training_loyo_skill"] = model.validation.global_skill
    return payload


__all__ = [
    "FINAL_TRAINING_PERIOD_LABEL",
    "RIDGE_LAMBDA_DEFAULT",
    "ErrorMetrics",
    "LeaveOneYearOutReport",
    "LinearModelFit",
    "ParameterLayout",
    "StochasticModelFit",
    "TrainingResult",
    "YearlyDesignStats",
    "YearlyValidationMetrics",
    "build_annual_basis",
    "build_global_linear_matrix",
    "build_parameter_payload",
    "compute_sufficient_stats",
    "fit_from_stats",
    "prepare_training_frame",
    "solve_normal_equations",
    "train_models",
    "train_stochastic_model",
]
