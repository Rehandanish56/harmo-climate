"""Core astronomical, thermodynamic, and linear-model helpers for HarmoClimate."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm, skewnorm

from .config import SAMPLES_PER_DAY
from .psychrometrics import (
    dew_point_c_from_e,
    relative_humidity_percent_from_specific,
    specific_humidity_kg_per_kg,
    thermo_from_T_P_RH,
    vapor_partial_pressure_hpa_from_q_p,
)
SOLAR_YEAR_DAYS: float = 365.242189
SOLAR_EPOCH_UTC = pd.Timestamp("2000-01-01 00:00:00", tz="UTC")
DATASET_COLUMNS: tuple[str, ...] = (
    "STATION_CODE",
    "STATION_NAME",
    "DT_UTC",
    "T",
    "RH",
    "P",
    "LON",
    "LAT",
    "ALTI",
)
_DERIVED_DATASET_COLUMNS: tuple[str, ...] = (
    "yday_frac_solar",
    "hour_solar",
    "delta_utc_solar_h",
    "Q",
    "Td",
    "E",
)
DEFAULT_PREPARED_COLUMNS: tuple[str, ...] = (
    "DT_UTC",
    "LON",
    "T",
    "RH",
    "P",
    *_DERIVED_DATASET_COLUMNS,
)
LINEAR_MODEL_VARIABLE_CODES: tuple[str, ...] = ("T", "RH", "TD", "Q", "E", "P")


def load_parquet_dataset(parquet_path: Path) -> pd.DataFrame:
    """Load a cached Parquet dataset and optionally constrain the returned columns."""

    if not parquet_path.exists():
        raise FileNotFoundError(f"No Parquet dataset found at {parquet_path}")

    df = pd.read_parquet(parquet_path)

    missing = set(DATASET_COLUMNS) - set(df.columns)
    if missing:
        raise KeyError(f"Dataset at {parquet_path} is missing required columns: {sorted(missing)}")

    return df


def compute_solar_time(
    dt_utc: pd.Series | np.ndarray,
    lon_deg: pd.Series | np.ndarray,
) -> pd.DataFrame:
    """
    Convert UTC datetimes and longitudes into solar descriptors.

    Returns:
        DataFrame with columns:
            - yday_frac_solar: solar-day index in a tropical year (no hour component),
            - hour_solar: local solar hour in [0, 24),
            - delta_utc_solar_h: UTC→solar offset in hours.
    """

    # Ensure inputs are Series for consistent handling
    if not isinstance(dt_utc, pd.Series):
        dt_utc = pd.Series(dt_utc)
    
    if not isinstance(lon_deg, pd.Series):
        lon_deg = pd.Series(lon_deg)

    # Align indices if possible, but prioritize length check
    if len(lon_deg) != len(dt_utc):
        raise ValueError("Longitude and datetime inputs must share the same length.")
    
    # Create a working copy to avoid side effects
    dt_series = pd.to_datetime(dt_utc, errors="coerce")
    lon_series = pd.to_numeric(lon_deg, errors="coerce")
    
    # If indices don't match, reset them to align by position (common issue with mixed inputs)
    if not dt_series.index.equals(lon_series.index):
         lon_series.index = dt_series.index

    # Ensure UTC tz-awareness
    if dt_series.dt.tz is None:
        dt_series = dt_series.dt.tz_localize("UTC")
    else:
        dt_series = dt_series.dt.tz_convert("UTC")

    # 1. Solar Day Fraction
    # Calculate days since epoch
    dt_utc_date = dt_series.dt.floor("D")
    delta_days = (dt_utc_date - SOLAR_EPOCH_UTC).dt.total_seconds() / 86400.0
    
    # Adjust for longitude (local solar time approximation)
    # 360 degrees = 1 day, so lon/360 is the fraction of a day adjustment
    solar_day_raw = delta_days + (lon_series / 360.0)
    solar_day = np.mod(solar_day_raw, SOLAR_YEAR_DAYS)

    # 2. Solar Hour
    # UTC decimal hour
    hour_utc = (
        dt_series.dt.hour
        + dt_series.dt.minute / 60.0
        + dt_series.dt.second / 3600.0
    )
    
    # Offset in hours (15 degrees per hour)
    delta_utc_solar_h = lon_series / 15.0
    hour_solar = (hour_utc + delta_utc_solar_h) % 24.0

    return pd.DataFrame(
        {
            "yday_frac_solar": solar_day.astype(np.float32),
            "hour_solar": hour_solar.astype(np.float32),
            "delta_utc_solar_h": delta_utc_solar_h.astype(np.float32),
        },
        index=dt_series.index,
    )


def prepare_dataset(
    df: pd.DataFrame,
    columns: Iterable[str] | None = None,
) -> pd.DataFrame:
    """
    Return a dataset copy populated with derived solar descriptors and moist thermodynamics.

    Args:
        df: Raw station dataset containing at least the core meteorological columns.
        columns: Optional iterable of column names to keep in the returned frame. Defaults to
            `DEFAULT_PREPARED_COLUMNS`.
    """

    required = {"DT_UTC", "LON", "T", "RH", "P"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Dataset is missing required columns: {sorted(missing)}")

    working = df.copy()
    working["DT_UTC"] = pd.to_datetime(working["DT_UTC"], utc=True, errors="coerce")
    working["LON"] = pd.to_numeric(working["LON"], errors="coerce")
    working = working.dropna(subset=["DT_UTC", "LON"])

    # Determine which columns to return
    if columns is None:
        requested = set(DEFAULT_PREPARED_COLUMNS)
    else:
        requested = set(columns)

    # Validate requested columns are either available or derivable
    available_base = set(working.columns)
    derived_names = set(_DERIVED_DATASET_COLUMNS)
    unknown = requested - (available_base | derived_names)
    if unknown:
        raise KeyError(f"Requested columns are not available: {sorted(unknown)}")

    # Identify what needs to be computed
    # Solar columns
    solar_cols = {"yday_frac_solar", "hour_solar", "delta_utc_solar_h"}
    if requested & solar_cols:
        solar_time = compute_solar_time(working["DT_UTC"], working["LON"])
        working = pd.concat([working, solar_time], axis=1)

    # Thermodynamic columns
    # Q depends on T, RH, P
    if "Q" in requested or "E" in requested or "Td" in requested:
         working["Q"] = specific_humidity_kg_per_kg(working["T"], working["RH"], working["P"])

    # E depends on Q, P (or T, RH, P indirectly)
    if "E" in requested or "Td" in requested:
         working["E"] = vapor_partial_pressure_hpa_from_q_p(working["Q"], working["P"])

    # Td depends on E
    if "Td" in requested:
         working["Td"] = dew_point_c_from_e(working["E"])

    # Final check to ensure all requested columns are present (e.g. if calculation failed silently)
    missing_requested = requested - set(working.columns)
    if missing_requested:
        raise KeyError(f"Unable to populate requested columns: {sorted(missing_requested)}")

    return working.loc[:, list(requested)].copy()


def load_stochastic_model(json_path: Path) -> dict:
    """Load a stochastic model JSON and validate required metadata."""
    with open(json_path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    if "metadata" not in payload or "model" not in payload:
        raise ValueError(f"{json_path.name} must contain 'metadata' and 'model' objects.")

    model = payload["model"]
    if "coefficient_matrix" not in model or "covariance_matrix" not in model:
        raise ValueError(
            f"{json_path.name} must include model.coefficient_matrix and model.covariance_matrix."
        )

    return payload


# ----------------------------- Linear model simulation -----------------------------

def _parameter_map(model_payload: dict) -> Dict[str, Tuple[np.ndarray, int]]:
    """Return a mapping of parameter name -> (coefficients, n_annual)."""

    model = model_payload["model"]
    coefficients = model.get("coefficients", [])
    layout = model.get("params_layout", [])

    mapping: Dict[str, Tuple[np.ndarray, int]] = {}
    for entry in layout:
        name = entry["name"]
        start = int(entry["start"])
        length = int(entry["length"])
        n_annual = int(entry["n_annual"])
        coeff_slice = np.asarray(coefficients[start : start + length], dtype=np.float64)
        mapping[name] = (coeff_slice, n_annual)
    return mapping


def _wrap_day(day: float) -> float:
    period = SOLAR_YEAR_DAYS
    day %= period
    if day < 0.0:
        day += period
    return day


def _wrap_hour(hour: float) -> float:
    hour %= 24.0
    if hour < 0.0:
        hour += 24.0
    return hour


def _eval_annual(coeffs: np.ndarray, n_annual: int, day: float) -> float:
    omega = 2.0 * math.pi / SOLAR_YEAR_DAYS
    value = float(coeffs[0])
    for k in range(1, n_annual + 1):
        angle = k * omega * day
        value += float(coeffs[2 * k - 1]) * math.cos(angle)
        value += float(coeffs[2 * k]) * math.sin(angle)
    return value


def _eval_annual_series(coeffs: Sequence[float], day_series: np.ndarray) -> np.ndarray:
    """Evaluate harmonic series for an array of days."""
    n_coeffs = len(coeffs)
    n_annual = (n_coeffs - 1) // 2
    omega = 2.0 * math.pi / SOLAR_YEAR_DAYS
    
    # Start with constant term
    values = np.full_like(day_series, coeffs[0], dtype=float)
    
    for k in range(1, n_annual + 1):
        angle = k * omega * day_series
        values += coeffs[2 * k - 1] * np.cos(angle)
        values += coeffs[2 * k] * np.sin(angle)
        
    return values


def skew_params_from_a_sigma(a: np.ndarray, sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Given shape parameter a and effective std sigma, compute (loc, scale)
    for a skew-normal distribution with mean 0 and std sigma.
    """
    a = np.asarray(a, dtype=float)
    sigma = np.asarray(sigma, dtype=float)

    delta = a / np.sqrt(1.0 + a * a)
    m1 = delta * math.sqrt(2.0 / math.pi)
    v1 = 1.0 - 2.0 * delta * delta / math.pi

    eps = 1e-12
    v1 = np.maximum(v1, eps)

    scale = sigma / np.sqrt(v1)
    loc = -scale * m1
    return loc, scale


def predict_model_solar(model_payload: dict, day_solar: float, hour_solar: float) -> float:
    """Evaluate a linear harmonic model in solar coordinates."""

    mapping = _parameter_map(model_payload)
    if "c0" not in mapping:
        raise ValueError("Model payload does not include mandatory parameter 'c0'.")

    day = math.floor(_wrap_day(day_solar))
    max_day_index = int(math.floor(SOLAR_YEAR_DAYS))
    if day > max_day_index:
        day = max_day_index
    hour = _wrap_hour(hour_solar)

    omega_hour = 2.0 * math.pi / 24.0
    result = _eval_annual(*mapping["c0"], day)

    n_diurnal = int(model_payload["model"].get("n_diurnal", 0))
    for m in range(1, n_diurnal + 1):
        cos_name = f"a{m}"
        sin_name = f"b{m}"
        if cos_name in mapping:
            result += _eval_annual(*mapping[cos_name], day) * math.cos(m * omega_hour * hour)
        if sin_name in mapping:
            result += _eval_annual(*mapping[sin_name], day) * math.sin(m * omega_hour * hour)
    return result


def climate_predict_solar(
    day_solar: float,
    hour_solar: float,
    temperature_model: dict,
    specific_humidity_model: dict,
    pressure_model: dict,
) -> Tuple[float, float, float]:
    """Return (temperature °C, specific humidity kg/kg, pressure hPa) in solar coordinates."""

    temp_c = predict_model_solar(temperature_model, day_solar, hour_solar)
    q = predict_model_solar(specific_humidity_model, day_solar, hour_solar)
    p_hpa = predict_model_solar(pressure_model, day_solar, hour_solar)
    return temp_c, q, p_hpa


def model_daily_stats_one_year_factorized(
    temperature_model: dict,
    specific_humidity_model: dict,
    pressure_model: dict,
    n_days: int = int(math.floor(SOLAR_YEAR_DAYS)),
    samples_per_day: int = SAMPLES_PER_DAY,
) -> Tuple[np.ndarray, Dict[str, Dict[str, np.ndarray]]]:
    """Compute model-driven daily envelopes for one full year."""

    hours = np.linspace(0.0, 24.0, samples_per_day, endpoint=False)
    days = np.arange(1, n_days + 1, dtype=int)

    stats: Dict[str, Dict[str, np.ndarray]] = {
        var: {
            "min": np.empty(n_days),
            "max": np.empty(n_days),
            "avg": np.empty(n_days),
        }
        for var in LINEAR_MODEL_VARIABLE_CODES
    }

    for idx, d in enumerate(days):
        temps: list[float] = []
        q_vals: list[float] = []
        pressures: list[float] = []
        for h in hours:
            temp_c, q, p_hpa = climate_predict_solar(
                float(d),
                float(h),
                temperature_model,
                specific_humidity_model,
                pressure_model,
            )
            temps.append(temp_c)
            q_vals.append(q)
            pressures.append(p_hpa)

        T_array = np.asarray(temps)
        Q_array = np.asarray(q_vals)
        P_array = np.asarray(pressures)
        RH_array = relative_humidity_percent_from_specific(T_array, Q_array, P_array)
        E_array = vapor_partial_pressure_hpa_from_q_p(Q_array, P_array)
        Td_array = dew_point_c_from_e(E_array)

        def store(code: str, arr: np.ndarray):
            stats[code]["min"][idx] = float(arr.min())
            stats[code]["max"][idx] = float(arr.max())
            stats[code]["avg"][idx] = float(arr.mean())

        store("T", T_array)
        store("Q", Q_array)
        store("P", P_array)
        store("RH", RH_array)
        store("E", E_array)
        store("TD", Td_array)

    return days, stats


def model_intraday_solar(
    temperature_model: dict,
    specific_humidity_model: dict,
    pressure_model: dict,
    day_solar: int,
    samples_per_day: int = SAMPLES_PER_DAY,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Return intraday solar-time profiles for the given models and day."""

    hours = np.linspace(0.0, 24.0, samples_per_day, endpoint=False)
    temps, qs, ps = [], [], []
    for h in hours:
        temp_c, q, p = climate_predict_solar(
            float(day_solar),
            float(h),
            temperature_model,
            specific_humidity_model,
            pressure_model,
        )
        temps.append(temp_c)
        qs.append(q)
        ps.append(p)

    temps_arr = np.asarray(temps, dtype=float)
    qs_arr = np.asarray(qs, dtype=float)
    ps_arr = np.asarray(ps, dtype=float)
    e_arr = vapor_partial_pressure_hpa_from_q_p(qs_arr, ps_arr)
    td_arr = dew_point_c_from_e(e_arr)
    rh_arr = relative_humidity_percent_from_specific(temps_arr, qs_arr, ps_arr)

    data = {
        "T": temps_arr,
        "Q": qs_arr,
        "P": ps_arr,
        "E": e_arr,
        "TD": td_arr,
        "RH": rh_arr,
    }
    return hours, data


def simulate_stochastic_year(
    temperature_model: dict,
    specific_humidity_model: dict,
    pressure_model: dict,
    stochastic_model: dict,
    year: int = 2010,
    samples_per_day: int = 24,
    seed: int | None = 42,
) -> pd.DataFrame:
    """Generate a one-year stochastic simulation combining harmonic baseline with VAR1 residuals.
    
    Args:
        temperature_model: Harmonic model for temperature
        specific_humidity_model: Harmonic model for specific humidity
        pressure_model: Harmonic model for pressure
        stochastic_model: VAR1 stochastic model with coefficient_matrix and covariance_matrix
        year: Year to simulate (default: 2025)
        samples_per_day: Samples per day, typically 24 for hourly data
        seed: Random seed for reproducibility (default: 42, None for random)
    
    Returns:
        DataFrame with simulated values and baseline for T, Q, P and derived variables
    """
    # Generate timestamps for the year
    start_date = f"{year}-01-01"
    end_date = f"{year+1}-01-01"
    timestamps = pd.date_range(start_date, end_date, freq="h", inclusive="left")
    
    # Compute solar time coordinates
    n_steps = len(timestamps)
    hours_total = np.arange(n_steps, dtype=float)
    day_solar = (hours_total / 24.0) % SOLAR_YEAR_DAYS + 1.0
    hour_solar = hours_total % 24.0
    
    # Compute harmonic baseline for T, Q, P
    T_base = np.zeros(n_steps)
    Q_base = np.zeros(n_steps)
    P_base = np.zeros(n_steps)
    
    for i in range(n_steps):
        T_base[i] = predict_model_solar(temperature_model, day_solar[i], hour_solar[i])
        Q_base[i] = predict_model_solar(specific_humidity_model, day_solar[i], hour_solar[i])
        P_base[i] = predict_model_solar(pressure_model, day_solar[i], hour_solar[i])
    
    # Generate VAR(p) residuals: X_t = A @ X_{lag} + noise
    # Variable order: P, Q, T (as defined in training.py)
    A = np.array(stochastic_model["model"]["coefficient_matrix"])
    Sigma = np.array(stochastic_model["model"]["covariance_matrix"])
    
    # Determine lags from metadata
    lags = stochastic_model["metadata"].get("lags")
    if lags is None:
        raise ValueError("Stochastic model metadata must include 'lags'.")
            
    max_lag = max(lags)

    rng = np.random.default_rng(seed=seed)
    noise = rng.multivariate_normal(mean=np.zeros(3), cov=Sigma, size=n_steps)
    
    residuals = np.zeros((n_steps, 3))
    # History buffer: [X_{t-1}, X_{t-2}, ..., X_{t-max_lag}]
    # history[k] corresponds to lag k+1
    history = np.zeros((max_lag, 3))
    
    for i in range(n_steps):
        # Construct lag vector for sparse lags: concatenate X_{t-k} for k in lags
        lag_components = []
        for lag in lags:
            # history[lag-1] is X_{t-lag}
            lag_components.append(history[lag-1])
            
        lag_vector = np.concatenate(lag_components)
        x_t = lag_vector @ A + noise[i]
        residuals[i, :] = x_t

        # Update history: shift old lags down, insert new x_t at top
        history[1:] = history[:-1]
        history[0] = x_t
    
    # Extract residuals by variable
    # Note: These are Z-scores (Gaussian domain)
    P_z = residuals[:, 0]
    Q_z = residuals[:, 1]
    T_z = residuals[:, 2]
    
    # Transform Z-scores back to original domain using Skew Normal parameters
    # x = skewnorm.ppf(norm.cdf(z), a, loc, scale)
    skew_params = stochastic_model["model"].get("skew_normal_params", {})
    
    def inverse_transform(z_scores, var_name):
        if var_name not in skew_params:
            # Fallback if no skew params (should not happen with new models)
            return z_scores
            
        params = skew_params[var_name]
        params = skew_params[var_name]
        
        # New harmonic envelope: expect "a" and "sigma" coefficients
        if "a" in params and "sigma" in params:
            coeffs_a = params["a"]
            coeffs_sigma = params["sigma"]
            
            a_ts = _eval_annual_series(coeffs_a, day_solar)
            sigma_ts = _eval_annual_series(coeffs_sigma, day_solar)
            
            loc_ts, scale_ts = skew_params_from_a_sigma(a_ts, sigma_ts)
            
            # Transform Z-score to uniform [0, 1]
            u = norm.cdf(z_scores)
            
            # Clip for numerical stability
            epsilon = 1e-6
            u = np.clip(u, epsilon, 1.0 - epsilon)
            
            # Transform uniform to original domain using Skew Normal PPF
            x = skewnorm.ppf(u, a_ts, loc=loc_ts, scale=scale_ts)
            return x
            
        else:
            # Fallback or error if format is unexpected
            # We strictly expect "a" and "sigma" now, but let's be safe
            # If we really want to drop backward compat, we can just raise or return z_scores
            return z_scores

    P_res = inverse_transform(P_z, "P")
    Q_res = inverse_transform(Q_z, "Q")
    T_res = inverse_transform(T_z, "T")
    
    # Combine baseline + residuals
    T_sim = T_base + T_res
    Q_sim = Q_base + Q_res
    P_sim = P_base + P_res
    
    # Compute derived variables
    E_sim = vapor_partial_pressure_hpa_from_q_p(Q_sim, P_sim)
    Td_sim = dew_point_c_from_e(E_sim)
    RH_sim = relative_humidity_percent_from_specific(T_sim, Q_sim, P_sim)
    
    return pd.DataFrame({
        "DT_UTC": timestamps,
        "day_solar": day_solar,
        "hour_solar": hour_solar,
        "T": T_sim,
        "Q": Q_sim,
        "P": P_sim,
        "E": E_sim,
        "TD": Td_sim,
        "RH": RH_sim,
        "T_base": T_base,
        "Q_base": Q_base,
        "P_base": P_base,
    })


__all__ = [
    "SOLAR_YEAR_DAYS",
    "SOLAR_EPOCH_UTC",
    "DATASET_COLUMNS",
    "specific_humidity_kg_per_kg",
    "compute_solar_time",
    "load_parquet_dataset",
    "load_stochastic_model",
    "predict_model_solar",
    "climate_predict_solar",
    "model_daily_stats_one_year_factorized",
    "model_intraday_solar",
    "prepare_dataset",
    "relative_humidity_percent_from_specific",
    "simulate_stochastic_year",
    "dew_point_c_from_e",
    "vapor_partial_pressure_hpa_from_q_p",
    "thermo_from_T_P_RH",
    "skew_params_from_a_sigma",
]
