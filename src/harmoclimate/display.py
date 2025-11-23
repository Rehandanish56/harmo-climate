"""Visualization helpers for HarmoClimate linear harmonic models."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Sequence, Tuple, Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from .config import SAMPLES_PER_DAY
from .core import (
    SOLAR_YEAR_DAYS,
    climate_predict_solar,
    load_parquet_dataset,
    load_stochastic_model,
    model_daily_stats_one_year_factorized,
    model_intraday_solar,
    prepare_dataset,
    predict_model_solar,
    simulate_stochastic_year,
)
from .psychrometrics import (
    dew_point_c_from_e,
    relative_humidity_percent_from_specific,
    vapor_partial_pressure_hpa_from_q_p,
)

MODEL_LINEWIDTH = 1.8
HISTORY_LINEWIDTH = 1.5
GRID_LINEWIDTH = 0.5
DISPLAY_VARIABLE_CHOICES = ("T", "RH", "TD", "Q", "E", "P")
DISPLAY_VARIABLE_DEFAULT = ("T", "Q", "P")


@dataclass(frozen=True)
class VariableConfig:
    title: str
    unit: str
    scale: float = 1.0


VARIABLE_METADATA = {
    "T": VariableConfig("Temperature", "°C"),
    "RH": VariableConfig("Relative humidity", "%"),
    "TD": VariableConfig("Dew point", "°C"),
    "Q": VariableConfig("Specific humidity", "g/kg", scale=1000.0),
    "E": VariableConfig("Vapor pressure", "hPa"),
    "P": VariableConfig("Pressure", "hPa"),
}


def normalize_display_variables(variables: Sequence[str] | None) -> Tuple[str, ...]:
    """Return a sanitized tuple of display variable codes in user-specified order."""

    if variables is None:
        return DISPLAY_VARIABLE_DEFAULT

    normalized: list[str] = []
    seen: set[str] = set()
    for entry in variables:
        code = str(entry).strip().upper()
        if not code:
            continue
        if code not in DISPLAY_VARIABLE_CHOICES:
            raise ValueError(
                f"Unsupported variable '{entry}'. Expected one of {', '.join(DISPLAY_VARIABLE_CHOICES)}."
            )
        if code not in seen:
            normalized.append(code)
            seen.add(code)

    if not normalized:
        return DISPLAY_VARIABLE_DEFAULT
    return tuple(normalized)


def load_linear_model(json_path: Path) -> dict:
    """Load a single-target linear model JSON and validate required metadata."""

    with open(json_path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    if "metadata" not in payload or "model" not in payload:
        raise ValueError(f"{json_path.name} must contain 'metadata' and 'model' objects.")

    metadata = payload["metadata"]
    required_meta = {
        "author",
        "version",
        "generated_at_utc",
        "country_code",
        "station_usual_name",
        "longitude_deg",
        "delta_utc_solar_h",
        "target_variable",
        "target_unit",
    }
    missing_meta = required_meta - metadata.keys()
    if missing_meta:
        raise ValueError(f"Missing required metadata fields: {sorted(missing_meta)}")

    model = payload["model"]
    if "coefficients" not in model or "params_layout" not in model:
        raise ValueError(f"{json_path.name} must include model.coefficients and model.params_layout.")

    return payload


def load_history_from_sample_data(parquet_path: Path) -> pd.DataFrame:
    """Load the cached dataset and emit key columns for comparison plots."""

    df = load_parquet_dataset(parquet_path)
    df = prepare_dataset(
        df,
        columns=(
            "yday_frac_solar",
            "hour_solar",
            "T",
            "RH",
            "Td",
            "P",
            "Q",
            "E",
            "DT_UTC",
        ),
    )
    df = df.dropna(
        subset=[
            "yday_frac_solar",
            "hour_solar",
            "T",
            "RH",
            "Td",
            "P",
            "Q",
            "E",
        ],
    )

    # Ensure day indices used in displays are 1-based to align with synthetic simulations
    # (simulate_stochastic_year starts at day_solar=1).
    max_day_index = int(math.floor(SOLAR_YEAR_DAYS))
    df["day"] = np.floor(df["yday_frac_solar"]).astype(int) + 1
    df.loc[df["day"] > max_day_index, "day"] = max_day_index
    df["hour"] = np.mod(df["hour_solar"].astype(float), 24.0)
    return df


def historical_climatology_daily(df: pd.DataFrame):
    """Aggregate historical daily envelopes for temperature, RH, Td, Q, E, and pressure."""

    grouped = df.groupby(["day", "hour"], as_index=False).agg(
        T_hour_mean=("T", "mean"),
        RH_hour_mean=("RH", "mean"),
        Td_hour_mean=("Td", "mean"),
        Q_hour_mean=("Q", "mean"),
        E_hour_mean=("E", "mean"),
        P_hour_mean=("P", "mean"),
    )
    daily = grouped.groupby("day", as_index=False).agg(
        temp_min_c=("T_hour_mean", "min"),
        temp_max_c=("T_hour_mean", "max"),
        temp_mean_c=("T_hour_mean", "mean"),
        rh_min_percent=("RH_hour_mean", "min"),
        rh_max_percent=("RH_hour_mean", "max"),
        rh_mean_percent=("RH_hour_mean", "mean"),
        td_min_c=("Td_hour_mean", "min"),
        td_max_c=("Td_hour_mean", "max"),
        td_mean_c=("Td_hour_mean", "mean"),
        q_min_kg_kg=("Q_hour_mean", "min"),
        q_max_kg_kg=("Q_hour_mean", "max"),
        q_mean_kg_kg=("Q_hour_mean", "mean"),
        e_min_hpa=("E_hour_mean", "min"),
        e_max_hpa=("E_hour_mean", "max"),
        e_mean_hpa=("E_hour_mean", "mean"),
        p_min_hpa=("P_hour_mean", "min"),
        p_max_hpa=("P_hour_mean", "max"),
        p_mean_hpa=("P_hour_mean", "mean"),
    )
    full = pd.DataFrame({"day": np.arange(1, int(math.floor(SOLAR_YEAR_DAYS)) + 1)})
    daily = full.merge(daily, on="day", how="left").interpolate(limit_direction="both")
    return daily, grouped


def _plot_variable_on_axis(
    ax: plt.Axes,
    code: str,
    x_data: np.ndarray,
    model_data: Dict[str, Any] | np.ndarray,
    hist_data: Dict[str, Any] | np.ndarray | None = None,
    x_data_hist: np.ndarray | None = None,
    is_intraday: bool = False,
):
    config = VARIABLE_METADATA[code]
    scale = config.scale
    
    ax.grid(True, which="both", linewidth=GRID_LINEWIDTH, alpha=0.4)
    ax.set_ylabel(f"({config.unit})")
    ax.set_title(config.title, loc="center", fontsize=11)

    if is_intraday:
        # model_data is array
        ax.plot(x_data, model_data * scale, label="Model", linewidth=MODEL_LINEWIDTH)
        if hist_data is not None and len(hist_data) > 0:
            x_hist = x_data_hist if x_data_hist is not None else x_data
            ax.plot(x_hist, hist_data * scale, "--", label="History mean", linewidth=HISTORY_LINEWIDTH, color="black")
    else:
        # model_data is dict with min/max/avg
        mn = model_data["min"] * scale
        mx = model_data["max"] * scale
        avg = model_data["avg"] * scale
        ax.fill_between(x_data, mn, mx, alpha=0.20, label="Model (min-max)")
        ax.plot(x_data, avg, linewidth=MODEL_LINEWIDTH, label="Model mean")
        
        if hist_data is not None:
             # hist_data is dict with min/max/mean
            h_mn = hist_data["min"] * scale
            h_mx = hist_data["max"] * scale
            h_avg = hist_data["mean"] * scale
            ax.fill_between(x_data, h_mn, h_mx, alpha=0.10, color="gray", label="History (min-max)")
            ax.plot(x_data, h_avg, linewidth=HISTORY_LINEWIDTH, linestyle="--", color="black", label="History mean")

    ax.legend(loc="best" if is_intraday else "upper right")


def plot_year(
    days: np.ndarray,
    model_stats: Dict[str, Dict[str, np.ndarray]],
    station_name: str,
    save_path: Path | str | None = None,
    hist_daily: pd.DataFrame | None = None,
    *,
    variables: Sequence[str] | None = None,
) -> Path | Figure:
    """Render annual envelopes for selected variables and optional history overlays."""

    enabled = normalize_display_variables(variables)
    axes_count = len(enabled)
    height = 4.0 + 1.8 * axes_count
    fig, axes = plt.subplots(axes_count, 1, figsize=(12, height), sharex=True)
    axes_list = [axes] if axes_count == 1 else list(axes)
    axis_map = {code: axes_list[idx] for idx, code in enumerate(enabled)}
    fig.suptitle(f"{station_name} — Annual model envelopes", fontsize=14)

    hist_stats = {}
    if hist_daily is not None and not hist_daily.empty:
        hist_sorted = hist_daily.sort_values("day")
        # Pre-process history into dict structure to match model_stats structure roughly
        # Mapping from our variable codes to dataframe columns
        col_map = {
            "T": "temp",
            "RH": "rh",
            "TD": "td",
            "Q": "q",
            "E": "e",
            "P": "p"
        }
        unit_suffix_map = {
            "T": "_c",
            "RH": "_percent",
            "TD": "_c",
            "Q": "_kg_kg",
            "E": "_hpa",
            "P": "_hpa"
        }
        
        for code in enabled:
            base = col_map.get(code)
            suffix = unit_suffix_map.get(code)
            if base and suffix:
                hist_stats[code] = {
                    "min": hist_sorted[f"{base}_min{suffix}"].to_numpy(dtype=float),
                    "max": hist_sorted[f"{base}_max{suffix}"].to_numpy(dtype=float),
                    "mean": hist_sorted[f"{base}_mean{suffix}"].to_numpy(dtype=float),
                }

    for code in enabled:
        _plot_variable_on_axis(
            axis_map[code],
            code,
            days,
            model_stats[code],
            hist_stats.get(code),
            x_data_hist=days,
            is_intraday=False
        )

    axes_list[-1].set_xlabel("Day of year")
    plt.tight_layout(rect=(0, 0, 1, 0.97))

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
        print(f"Figure saved -> {save_path}")
        return save_path

    plt.show()
    return fig


def history_intraday_mean(hourly_hist: pd.DataFrame, day: int):
    sub = hourly_hist[hourly_hist["day"] == day].sort_values("hour")
    if sub.empty:
        return (
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
        )
    return (
        sub["hour"].to_numpy(),
        sub["T_hour_mean"].to_numpy(),
        sub["RH_hour_mean"].to_numpy(),
        sub["Td_hour_mean"].to_numpy(),
        sub["Q_hour_mean"].to_numpy(),
        sub["E_hour_mean"].to_numpy(),
        sub["P_hour_mean"].to_numpy(),
    )


def plot_intraday(
    temperature_model: dict,
    specific_humidity_model: dict,
    pressure_model: dict,
    hourly_hist: pd.DataFrame | None,
    station_name: str,
    samples_per_day: int = SAMPLES_PER_DAY,
    *,
    day: int = 15,
    save_path: Path | str | None = None,
    variables: Sequence[str] | None = None,
):
    """Render a static intraday plot for a specific solar day.

    When ``save_path`` is provided, the figure is written to disk and the resulting
    ``Path`` is returned; otherwise a Matplotlib figure is shown and returned.
    """

    max_day_index = int(math.floor(SOLAR_YEAR_DAYS))
    if day < 1 or day > max_day_index:
        raise ValueError(f"day must be within [1, {max_day_index}] (received {day}).")

    h_mod, model_data = model_intraday_solar(
        temperature_model,
        specific_humidity_model,
        pressure_model,
        day,
        samples_per_day=samples_per_day,
    )

    hist_data = {}
    h_hist = np.array([])
    if hourly_hist is not None:
        # history_intraday_mean returns tuple: (h, T, RH, Td, Q, E, P)
        # We need to map this to our dict structure
        h_hist, T_hist, RH_hist, Td_hist, Q_hist, E_hist, P_hist = history_intraday_mean(hourly_hist, day)
        if len(h_hist) > 0:
            hist_data = {
                "T": T_hist,
                "RH": RH_hist,
                "TD": Td_hist,
                "Q": Q_hist,
                "E": E_hist,
                "P": P_hist,
            }

    enabled = normalize_display_variables(variables)
    axes_count = len(enabled)
    height = 4.0 + 1.5 * axes_count
    fig, axes = plt.subplots(axes_count, 1, figsize=(9, height), sharex=True)
    axes_list = [axes] if axes_count == 1 else list(axes)
    axis_map = {code: axes_list[idx] for idx, code in enumerate(enabled)}
    fig.suptitle(f"{station_name} — Intraday profile (solar time) — day {day}", fontsize=12)

    for code in enabled:
        ax = axis_map[code]
        ax.set_xlim(0, 24)
        _plot_variable_on_axis(
            ax,
            code,
            h_mod,
            model_data[code],
            hist_data.get(code),
            x_data_hist=h_hist,
            is_intraday=True
        )

    axes_list[-1].set_xlabel("Solar hour")
    plt.tight_layout()

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
        print(f"Figure saved -> {save_path}")
        return save_path

    plt.show()
    return fig


def plot_stochastic_year(
    sim_df: pd.DataFrame,
    station_name: str,
    save_path: Path | str | None = None,
    real_df: pd.DataFrame | None = None,
    variables: Sequence[str] | None = None,
) -> Path | Figure:
    """Render stochastic simulation vs real data for a full year."""
    
    enabled = normalize_display_variables(variables)
    axes_count = len(enabled)
    height = 4.0 + 1.8 * axes_count
    fig, axes = plt.subplots(axes_count, 1, figsize=(12, height), sharex=True)
    axes_list = [axes] if axes_count == 1 else list(axes)
    axis_map = {code: axes_list[idx] for idx, code in enumerate(enabled)}
    fig.suptitle(f"{station_name} — Stochastic Simulation vs Real Data", fontsize=14)
    
    # X-axis: day of year
    sim_doy = sim_df["day_solar"]
    
    real_doy = None
    sim_year = int(sim_df["DT_UTC"].dt.year.iloc[0]) if "DT_UTC" in sim_df.columns else None
    real_label = f"Real Data ({sim_year})" if sim_year is not None else "Real Data"

    if real_df is not None:
        # Prefer aligning real data to the simulation year; fall back to latest year if absent
        if "DT_UTC" in real_df.columns:
            years = real_df["DT_UTC"].dt.year.unique()
            target_year = sim_year if sim_year in years else years.max()
            mask = real_df["DT_UTC"].dt.year == target_year
            real_df_year = real_df[mask]
            if not real_df_year.empty:
                real_df = real_df_year
                if sim_year != target_year:
                    real_label = f"Real Data ({target_year})"

        real_df = real_df.copy()
        if "yday_frac_solar" in real_df.columns:
            real_doy = real_df["yday_frac_solar"] + 1.0
        elif "DT_UTC" in real_df.columns:
            real_doy = real_df["DT_UTC"].dt.dayofyear + real_df["DT_UTC"].dt.hour / 24.0
    
    for code in enabled:
        ax = axis_map[code]
        config = VARIABLE_METADATA[code]
        scale = config.scale
        
        ax.grid(True, which="both", linewidth=GRID_LINEWIDTH, alpha=0.4)
        ax.set_ylabel(f"({config.unit})")
        ax.set_title(config.title, loc="center", fontsize=11)
        
        # Plot Simulation
        ax.plot(sim_doy, sim_df[code] * scale, label="Simulated", alpha=0.7, linewidth=0.8)
        
        # Plot Baseline (if available in sim_df, we added T_base etc)
        base_col = f"{code}_base"
        if base_col in sim_df.columns:
             ax.plot(sim_doy, sim_df[base_col] * scale, label="Harmonic Baseline", color="k", linestyle="--", alpha=0.6, linewidth=1)
        elif code in ["E", "TD", "RH"]:
             # These are derived, so we don't have direct baseline columns in sim_df easily unless we computed them.
             # For now, skip baseline for derived vars or compute them if needed.
             pass

        # Plot Real Data at native granularity to match simulation
        if real_df is not None and real_doy is not None:
            # Map code to column name in real_df
            # real_df usually has T, Q, P, RH, E, Td
            col_name = code
            if code == "TD" and "Td" in real_df.columns:
                col_name = "Td"
             
            if col_name in real_df.columns:
                ax.plot(real_doy, real_df[col_name] * scale, label=real_label, alpha=0.4, linewidth=0.8)
        
        ax.legend(loc="upper right")
        
    axes_list[-1].set_xlabel("Day of year")
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
        print(f"Figure saved -> {save_path}")
        return save_path

    plt.show()
    return fig


def plot_stochastic_intraday(
    sim_df: pd.DataFrame,
    station_name: str,
    day: int,
    save_path: Path | str | None = None,
    real_df: pd.DataFrame | None = None,
    variables: Sequence[str] | None = None,
) -> Path | Figure:
    """Render stochastic simulation vs real data for a specific day."""
    
    # Filter for the specific day; simulated day_solar is 1-based.
    sim_day = sim_df[sim_df["day_solar"].astype(int) == day]
    
    real_day = None
    if real_df is not None:
        # Align real data year to simulation year when possible; otherwise use latest.
        if "DT_UTC" in real_df.columns:
            years = real_df["DT_UTC"].dt.year.unique()
            target_year = sim_df["DT_UTC"].dt.year.iloc[0] if "DT_UTC" in sim_df.columns else years.max()
            if target_year not in years:
                target_year = years.max()
            real_df = real_df[real_df["DT_UTC"].dt.year == target_year]

        # Filter real data for the same day (approximate if years differ)
        if "yday_frac_solar" in real_df.columns:
             real_day = real_df[(np.floor(real_df["yday_frac_solar"]).astype(int) + 1) == day]
        elif "DT_UTC" in real_df.columns:
             real_day = real_df[real_df["DT_UTC"].dt.dayofyear == day]
    
    enabled = normalize_display_variables(variables)
    axes_count = len(enabled)
    height = 4.0 + 1.5 * axes_count
    fig, axes = plt.subplots(axes_count, 1, figsize=(9, height), sharex=True)
    axes_list = [axes] if axes_count == 1 else list(axes)
    axis_map = {code: axes_list[idx] for idx, code in enumerate(enabled)}
    fig.suptitle(f"{station_name} — Stochastic Intraday (Day {day})", fontsize=12)
    
    for code in enabled:
        ax = axis_map[code]
        config = VARIABLE_METADATA[code]
        scale = config.scale
        
        ax.grid(True, which="both", linewidth=GRID_LINEWIDTH, alpha=0.4)
        ax.set_ylabel(f"({config.unit})")
        ax.set_title(config.title, loc="center", fontsize=11)
        ax.set_xlim(0, 24)
        
        # Plot Simulation
        ax.plot(sim_day["hour_solar"], sim_day[code] * scale, label="Simulated", linewidth=1.5)
        
        # Plot Baseline
        base_col = f"{code}_base"
        if base_col in sim_day.columns:
             ax.plot(sim_day["hour_solar"], sim_day[base_col] * scale, label="Baseline", color="k", linestyle="--", alpha=0.6)
        
        # Plot Real Data
        if real_day is not None and not real_day.empty:
            col_name = code
            if code == "TD" and "Td" in real_day.columns:
                col_name = "Td"
            
            x_real = None
            if "hour_solar" in real_day.columns:
                x_real = real_day["hour_solar"]
            elif "DT_UTC" in real_day.columns:
                x_real = real_day["DT_UTC"].dt.hour + real_day["DT_UTC"].dt.minute / 60.0
            
            if col_name in real_day.columns and x_real is not None:
                ax.plot(x_real, real_day[col_name] * scale, label="Real", alpha=0.6, linewidth=1.5)
        
        ax.legend(loc="best")
        
    axes_list[-1].set_xlabel("Solar hour")
    plt.tight_layout()
    
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
        print(f"Figure saved -> {save_path}")
        return save_path

    plt.show()
    return fig


__all__ = [
    "DISPLAY_VARIABLE_CHOICES",
    "DISPLAY_VARIABLE_DEFAULT",
    "climate_predict_solar",
    "historical_climatology_daily",
    "history_intraday_mean",
    "load_linear_model",
    "load_history_from_sample_data",
    "model_daily_stats_one_year_factorized",
    "model_intraday_solar",
    "normalize_display_variables",
    "plot_intraday",
    "plot_year",
    "plot_stochastic_year",
    "plot_stochastic_intraday",
    "predict_model_solar",
]
