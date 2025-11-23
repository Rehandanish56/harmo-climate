"""End-to-end pipeline for training and exporting HarmoClimate models."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from .config import (
    ANNUAL_HARMONICS_PER_PARAM,
    AUTHOR_NAME,
    COUNTRY_CODE,
    DATA_DIR,
    DEFAULT_ANNUAL_HARMONICS,
    GENERATED_DIR,
    MEDIA_DIR,
    MODEL_DIR,
    N_DIURNAL_HARMONICS,
    STATION_CODE,
    TEMPLATE_DIR,
    ArtifactPaths,
    build_artifact_paths,
    build_urls_for_station,
    compute_output_basename,
    slugify_station_name,
)
from .core import (
    load_parquet_dataset,
    model_daily_stats_one_year_factorized,
    simulate_stochastic_year,
)
from .data_ingest import StationRecord, StreamResult, stream_filter_to_disk
from .metadata import StationMetadata, summarize_station
from .template_cpp import generate_cpp_header
from .evaluation import evaluate_loyo
from .training import (
    LinearModelFit,
    RIDGE_LAMBDA_DEFAULT,
    StochasticModelFit,
    build_parameter_payload,
    compute_sufficient_stats,
    prepare_training_frame,
    train_models,
    train_stochastic_model,
)
from .display import (
    historical_climatology_daily,
    load_history_from_sample_data,
    load_linear_model,
    load_stochastic_model,
    normalize_display_variables,
    plot_intraday,
    plot_year,
    plot_stochastic_year,
    plot_stochastic_intraday,
)
_TARGET_SUFFIXES = {
    "T": "_temperature",
    "Q": "_specific_humidity",
    "P": "_pressure",
}


def ensure_directories() -> None:
    for path in (GENERATED_DIR, DATA_DIR, MODEL_DIR, TEMPLATE_DIR, MEDIA_DIR):
        path.mkdir(parents=True, exist_ok=True)


def _export_training_metrics_files(model: LinearModelFit, model_path: Path) -> None:
    """Persist LOYO metrics alongside the primary model artefact."""

    report = getattr(model, "validation", None)
    if report is None:
        return

    metrics_dir = model_path.parent / "training_metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)

    metrics_stem = f"{model_path.stem}_training_metrics"
    metrics_json_path = metrics_dir / f"{metrics_stem}.json"
    metrics_csv_path = metrics_dir / f"{metrics_stem}.csv"

    payload = report.to_payload()

    with open(metrics_json_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)

    fieldnames = ["year", "rmse", "mse_model", "mse_ref", "skill", "n"]
    with open(metrics_csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        if report.years:
            for entry in report.years:
                row = entry.as_dict()
                writer.writerow(
                    {
                        "year": row["year"],
                        "rmse": row["rmse"],
                        "mse_model": row["mse_model"],
                        "mse_ref": row["mse_ref"],
                        "skill": row["skill"],
                        "n": row["n"],
                    }
                )

            total_n = sum(entry.n for entry in report.years)
            if total_n > 0:
                weighted_mse_model = sum(entry.n * entry.mse_model for entry in report.years)
                weighted_mse_ref = sum(entry.n * entry.mse_ref for entry in report.years)
                global_rmse = payload["global"]["rmse"]
                global_skill = payload["global"]["skill"]
                global_mse_model = weighted_mse_model / total_n
                global_mse_ref = weighted_mse_ref / total_n
                writer.writerow(
                    {
                        "year": "global",
                        "rmse": global_rmse,
                        "mse_model": global_mse_model,
                        "mse_ref": global_mse_ref,
                        "skill": global_skill,
                        "n": total_n,
                    }
                )

    if report.years:
        print(f"[OK] LOYO metrics exported to {metrics_json_path} and {metrics_csv_path}")
    else:
        print(
            "[Info] LOYO validation produced no per-year entries; "
            f"exported summary scaffold to {metrics_json_path} and {metrics_csv_path}",
        )


def _export_stochastic_model(
    model: StochasticModelFit,
    base_path: Path,
    metadata: dict[str, object],
    generation_date: str,
) -> None:
    """Export the stochastic model to a JSON file."""
    
    # Construct path: {station_slug}_stochastic.json
    # base_path is like .../models/fr_bourges_temperature.json
    # We want .../models/fr_bourges_stochastic.json
    
    # We can derive the station slug from the base path or pass it in.
    # Using _slug_from_model_path logic might be safer if available, 
    # but here we can just use the parent directory and a constructed filename 
    # if we assume base_path follows the standard naming.
    # However, let's try to be robust.
    
    # Actually, let's just use the directory of base_path and the station code/slug from metadata if possible,
    # or derive from base_path stem.
    
    # Let's assume the caller handles the path generation or we do it here.
    # In _finalize_pipeline we have artifact_paths.
    
    # Let's define the payload first.
    payload = {
        "metadata": {
            "version": "1.0",
            "model_type": model.model_type,
            "author": metadata.get("author", "HarmoClimate"),
            "generated_at_utc": generation_date,
            "station_code": metadata.get("station_code"),
            "station_usual_name": metadata.get("station_usual_name"),
            "variables": model.variables,
            "lags": model.lags,
        },
        "model": {
            "coefficient_matrix": model.coefficient_matrix_list(),
            "covariance_matrix": model.covariance_matrix_list(),
            "skew_normal_params": model.skew_normal_params,
        },
    }
    
    with open(base_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    
    print(f"[OK] Stochastic model exported to {base_path}")


def clean_pipeline() -> list[Path]:
    """Remove cached Parquet datasets under the generated data directory."""

    if not DATA_DIR.exists():
        print(f"[Info] Data directory {DATA_DIR} does not exist; nothing to clean.")
        return []

    removed: list[Path] = []
    for parquet_path in DATA_DIR.rglob("*.parquet"):
        if parquet_path.is_file():
            parquet_path.unlink()
            removed.append(parquet_path)
            print(f"[OK] Removed {parquet_path}")

    if not removed:
        print(f"[Info] No Parquet files found under {DATA_DIR}.")

    return removed


def _iso_utc(value) -> str | None:
    """Convert a timestamp to an ISO 8601 string in UTC."""
    if value is None:
        return None
    if hasattr(value, "to_pydatetime"):
        value = value.to_pydatetime()
    if isinstance(value, datetime):
        if value.tzinfo is None:
            value = value.replace(tzinfo=timezone.utc)
        else:
            value = value.astimezone(timezone.utc)
        return value.isoformat()
    return None


def _finalize_pipeline(
    *,
    df: pd.DataFrame,
    artifact_paths: ArtifactPaths,
    station_records: Iterable[StationRecord],
    station_name: str,
    station_code: str | None,
) -> StationMetadata:
    """Train models, export artefacts, and summarise metadata."""

    station_meta = summarize_station(station_records, df, fallback_name=station_name)
    if station_meta.station_code is None and station_code:
        station_meta = replace(station_meta, station_code=str(station_code))

    print(f"[Info] Mean station longitude = {station_meta.longitude_deg:.5f} °")
    if not math.isnan(station_meta.latitude_deg):
        print(f"[Info] Mean station latitude  = {station_meta.latitude_deg:.5f} °")
    else:
        print("[Info] Mean station latitude  = (not available in source)")
    print(f"[Info] UTC -> solar offset     = {station_meta.delta_utc_solar_h:.5f} h")
    print(f"[Info] Station usual name      = {station_meta.station_name}")
    if station_meta.station_code:
        print(f"[Info] Station code           = {station_meta.station_code}")

    prepared = prepare_training_frame(df)
    ridge_lambda = RIDGE_LAMBDA_DEFAULT
    annual_overrides = dict(ANNUAL_HARMONICS_PER_PARAM)

    result = train_models(
        prepared,
        n_diurnal=N_DIURNAL_HARMONICS,
        default_n_annual=DEFAULT_ANNUAL_HARMONICS,
        annual_per_param=annual_overrides,
        ridge_lambda=ridge_lambda,
    )

    reference_template = {
        "type": "climatology_mean",
        "time_basis": "UTC",
        "calendar": "no-leap",
        "grouping": "utc_day_of_year × utc_hour",
        "exclusion": "held-out year",
        "hours_per_day": 24,
        "days_per_year": 365,
    }
    model_spec = {
        "n_diurnal": int(N_DIURNAL_HARMONICS),
        "default_n_annual": int(DEFAULT_ANNUAL_HARMONICS),
        "annual_per_param": {k: int(v) for k, v in annual_overrides.items()},
        "ridge_lambda": float(ridge_lambda),
    }

    evaluation_meta = {
        "evaluation_time_base": "UTC",
        "model_time_base": "solar",
        "baseline": "climatology_mean per (utc_day, utc_hour), LOYO",
    }

    model_fits = {
        "T": result.temperature_model,
        "Q": result.specific_humidity_model,
        "P": result.pressure_model,
    }

    for target, model in model_fits.items():
        stats = compute_sufficient_stats(
            prepared,
            target=target,
            n_diurnal=N_DIURNAL_HARMONICS,
            default_n_annual=DEFAULT_ANNUAL_HARMONICS,
            annual_per_param=annual_overrides,
        )
        report = evaluate_loyo(
            stats,
            ridge_lambda=ridge_lambda,
            reference_spec=dict(reference_template),
        )
        report.hyperparameters = {
            "model": dict(model_spec),
            "reference": dict(reference_template),
            **evaluation_meta,
        }
        model.validation = report

    temp_metrics = result.temperature_model.metrics
    q_metrics = result.specific_humidity_model.metrics
    p_metrics = result.pressure_model.metrics

    print(f"Global MAE Temperature       : {temp_metrics.mae:.2f} °C")
    print(f"Global MAE Specific Humidity : {q_metrics.mae:.4f} kg/kg")
    print(f"Global MAE Pressure          : {p_metrics.mae:.2f} hPa")

    source_data_utc_start = None
    source_data_utc_end = None
    if "DT_UTC" in df.columns:
        dt_series = df["DT_UTC"].dropna()

        if not dt_series.empty:
            source_data_utc_start = _iso_utc(dt_series.min())
            source_data_utc_end = _iso_utc(dt_series.max())

    metadata_payload = {
        "station_code": station_meta.station_code,
        "station_usual_name": station_meta.station_name,
        "longitude_deg": station_meta.longitude_deg,
        "latitude_deg": station_meta.latitude_deg,
        "altitude_m": station_meta.altitude_m,
        "delta_utc_solar_h": station_meta.delta_utc_solar_h,
        "source_data_utc_start": source_data_utc_start,
        "source_data_utc_end": source_data_utc_end,
    }

    generation_date = datetime.now(timezone.utc).isoformat()

    payloads: dict[str, dict] = {}
    artifact_paths.model_temperature_json.parent.mkdir(parents=True, exist_ok=True)
    for target, model in model_fits.items():
        payload = build_parameter_payload(model, metadata_payload, generation_date)
        payloads[target] = payload

    path_map = {
        "T": artifact_paths.model_temperature_json,
        "Q": artifact_paths.model_specific_humidity_json,
        "P": artifact_paths.model_pressure_json,
    }

    for target, payload in payloads.items():
        path = path_map[target]
        _export_training_metrics_files(model_fits[target], path)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
        print(f"[OK] {payload['metadata']['target_variable']} parameters exported to {path}")

    # Train and export stochastic model
    stochastic_payload = None
    print("[Info] Training stochastic residual model...")
    try:
        stochastic_model = train_stochastic_model(result, prepared)
        
        # Determine path for stochastic model
        # We use the same directory as other models
        model_dir = artifact_paths.model_temperature_json.parent
        
        # Derive base name from temperature model path to ensure consistency (e.g. fr_bourges)
        temp_stem = artifact_paths.model_temperature_json.stem
        temp_suffix = _TARGET_SUFFIXES["T"]
        if temp_stem.endswith(temp_suffix):
            base_name = temp_stem[:-len(temp_suffix)]
        else:
            base_name = temp_stem
            
        stochastic_json_path = model_dir / f"{base_name}_stochastic.json"
        
        stochastic_payload = {
            "metadata": {
                "version": "1.0",
                "model_type": stochastic_model.model_type,
                "author": metadata_payload.get("author", AUTHOR_NAME),
                "generated_at_utc": generation_date,
                "station_code": metadata_payload.get("station_code"),
                "station_usual_name": metadata_payload.get("station_usual_name"),
                "variables": stochastic_model.variables,
                "lags": stochastic_model.lags,
            },
            "model": {
                "coefficient_matrix": stochastic_model.coefficient_matrix_list(),
                "covariance_matrix": stochastic_model.covariance_matrix_list(),
                "skew_normal_params": stochastic_model.skew_normal_params,
            },
        }

        with open(stochastic_json_path, "w", encoding="utf-8") as handle:
            json.dump(stochastic_payload, handle, indent=2)
        print(f"[OK] Stochastic model exported to {stochastic_json_path}")

    except Exception as e:
        print(f"[Warn] Failed to train/export stochastic model: {e}")

    generate_cpp_header(
        payloads["T"],
        payloads["Q"],
        payloads["P"],
        artifact_paths.cpp_header,
        stochastic_payload=stochastic_payload,
    )

    return station_meta


def _as_float(value: object) -> float:
    """Best-effort float coercion."""

    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _record_from_metadata(metadata: dict[str, object]) -> StationRecord:
    """Build a synthetic station record from cached metadata."""

    station_code_raw = metadata.get("station_code")
    station_code = str(station_code_raw).strip() if station_code_raw is not None else None
    if station_code == "":
        station_code = None

    station_name = str(
        metadata.get("station_usual_name")
        or metadata.get("station_name")
        or metadata.get("station_usual_name".upper())
        or "",
    ).strip()

    return StationRecord(
        station_code=station_code,
        station_name=station_name,
        lon=_as_float(metadata.get("longitude_deg")),
        lat=_as_float(metadata.get("latitude_deg")),
        alti=_as_float(metadata.get("altitude_m")),
        delta_utc_solar_h=_as_float(metadata.get("delta_utc_solar_h")),
    )


def _resolve_model_path(model_json: str | Path) -> Path:
    """Resolve a model JSON path, falling back to the generated models directory."""

    candidate = Path(model_json)
    if candidate.exists():
        return candidate

    fallback = MODEL_DIR / candidate.name
    if fallback.exists():
        return fallback

    raise FileNotFoundError(f"Could not locate model JSON at {candidate} or {fallback}")


def _slug_from_model_path(model_path: Path) -> str:
    """Extract the station slug from a model file name."""

    basename = model_path.stem
    expected_prefix = f"{COUNTRY_CODE.lower()}_"
    if basename.startswith(expected_prefix):
        basename = basename[len(expected_prefix) :]

    for suffix in _TARGET_SUFFIXES.values():
        if basename.endswith(suffix):
            basename = basename[: -len(suffix)]
            break
    
    if basename.endswith("_stochastic"):
        basename = basename[:-len("_stochastic")]

    return basename


def _apply_model_path_override(
    base_paths: ArtifactPaths,
    model_path: Path,
    target_variable: str,
) -> ArtifactPaths:
    """Return artefact paths overriding model JSON locations when needed."""

    suffix = _TARGET_SUFFIXES.get(target_variable)
    detected_target = target_variable
    if suffix is None:
        for candidate_target, candidate_suffix in _TARGET_SUFFIXES.items():
            if model_path.stem.endswith(candidate_suffix):
                detected_target = candidate_target
                suffix = candidate_suffix
                break
        else:
            suffix = ""

    stem = model_path.stem
    base_stem = stem[:-len(suffix)] if suffix and stem.endswith(suffix) else stem
    if base_stem.endswith("_stochastic"):
        base_stem = base_stem[:-len("_stochastic")]

    model_dir = model_path.parent
    derived_paths = {
        target: model_dir / f"{base_stem}{suffix_value}{model_path.suffix}"
        for target, suffix_value in _TARGET_SUFFIXES.items()
    }
    
    if not model_path.stem.endswith("_stochastic"):
        derived_paths[detected_target] = model_path
    
    # If target is stochastic, we don't override a specific model path in the artifact paths
    # but we ensure the derived paths are correct for T, Q, P.
    # The ArtifactPaths struct doesn't have a stochastic path field, but that's fine.
    
    if (
        derived_paths["T"] == base_paths.model_temperature_json
        and derived_paths["Q"] == base_paths.model_specific_humidity_json
        and derived_paths["P"] == base_paths.model_pressure_json
    ):
        return base_paths

    return ArtifactPaths(
        parquet=base_paths.parquet,
        model_temperature_json=derived_paths["T"],
        model_specific_humidity_json=derived_paths["Q"],
        model_pressure_json=derived_paths["P"],
        cpp_header=base_paths.cpp_header,
    )


def _normalize_model_basename(candidate: str) -> str:
    """Return the canonical model basename (e.g. fr_bourges) from user input."""

    name = Path(candidate).name
    if name.endswith(".json"):
        name = name[:-5]

    for suffix in _TARGET_SUFFIXES.values():
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    
    if name.endswith("_stochastic"):
        name = name[:-len("_stochastic")]

    prefix = f"{COUNTRY_CODE.lower()}_"
    if name.startswith(prefix):
        return name

    slug = slugify_station_name(name)
    if not slug:
        raise ValueError(f"Unable to derive station slug from model name '{candidate}'.")
    return f"{prefix}{slug}"


def run_pipeline(station_code: str = STATION_CODE, urls: Sequence[str] | None = None) -> StationMetadata:
    """Execute the full training, export, and template generation workflow."""

    ensure_directories()

    resolved_urls = list(urls) if urls is not None else build_urls_for_station(station_code)
    stream_result: StreamResult = stream_filter_to_disk(resolved_urls, station_code=station_code)
    artifact_paths = build_artifact_paths(stream_result.station_slug)
    df = load_parquet_dataset(artifact_paths.parquet)

    return _finalize_pipeline(
        df=df,
        artifact_paths=artifact_paths,
        station_records=stream_result.station_records,
        station_name=stream_result.station_name,
        station_code=station_code,
    )


def generate_pipeline(station_code: str) -> StationMetadata:
    """Fetch remote data for a station code and run the end-to-end pipeline."""

    return run_pipeline(station_code=station_code)


def regenerate_pipeline(model_json: str | Path) -> StationMetadata:
    """Rebuild artefacts from a cached dataset based on an existing model JSON."""

    ensure_directories()

    model_path = _resolve_model_path(model_json)
    with open(model_path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    metadata = payload.get("metadata", {})
    target_variable_raw = str(metadata.get("target_variable") or "").strip()
    station_code_raw = metadata.get("station_code")
    station_code = str(station_code_raw).strip() if station_code_raw is not None else None
    if station_code == "":
        station_code = None

    station_name = str(metadata.get("station_usual_name") or station_code or model_path.stem).strip()
    station_slug_hint = slugify_station_name(station_name) if station_name else ""
    if not station_slug_hint:
        station_slug_hint = _slug_from_model_path(model_path)

    artifact_paths = build_artifact_paths(station_slug_hint)

    if target_variable_raw in _TARGET_SUFFIXES:
        target_variable = target_variable_raw
    else:
        detected = None
        for candidate_target, candidate_suffix in _TARGET_SUFFIXES.items():
            if model_path.stem.endswith(candidate_suffix):
                detected = candidate_target
                break
        target_variable = detected or "T"

    artifact_paths = _apply_model_path_override(artifact_paths, model_path, target_variable)

    station_records: list[StationRecord] = []

    if artifact_paths.parquet.exists():
        print(f"[Cache] Using cached dataset at {artifact_paths.parquet}")
        df = load_parquet_dataset(artifact_paths.parquet)
        if metadata:
            station_records.append(_record_from_metadata(metadata))
    else:
        if not station_code:
            raise RuntimeError(
                "The provided model JSON does not include a station_code; cannot refresh remote data.",
            )
        print("[Info] Cached dataset not found; streaming source archives again.")
        stream_result: StreamResult = stream_filter_to_disk(
            build_urls_for_station(station_code),
            station_code=station_code,
        )
        artifact_paths = build_artifact_paths(stream_result.station_slug)
        artifact_paths = _apply_model_path_override(artifact_paths, model_path, target_variable)
        df = load_parquet_dataset(artifact_paths.parquet)
        station_records.extend(stream_result.station_records)
        station_name = stream_result.station_name

    return _finalize_pipeline(
        df=df,
        artifact_paths=artifact_paths,
        station_records=station_records,
        station_name=station_name,
        station_code=station_code,
    )


def display_pipeline(
    model_json: str | Path,
    *,
    mode: str = "annual",
    day: int | None = None,
    variables: Sequence[str] | None = None,
) -> Path:
    """Render plots for a target bundle with the requested variable panels."""

    ensure_directories()

    mode_normalized = (mode or "annual").strip().lower()
    if mode_normalized not in {"annual", "intraday", "stochastic_annual", "stochastic_intraday"}:
        raise ValueError(f"Unsupported display mode '{mode}'. Expected 'annual', 'intraday', 'stochastic_annual', or 'stochastic_intraday'.")
    if "intraday" in mode_normalized:
        if day is None:
            raise ValueError(f"The '{mode_normalized}' mode requires a --day argument.")
        day_int = int(day)
    else:
        day_int = None

    model_path = _resolve_model_path(model_json)
    
    # Check if this is a stochastic model file
    is_stochastic_file = model_path.name.endswith("_stochastic.json")
    
    if is_stochastic_file:
        # If user passed a stochastic file, we assume they want stochastic display unless specified otherwise
        # But if they passed mode="annual", we might want to switch to "stochastic_annual" automatically?
        # Or just proceed.
        if mode_normalized == "annual":
            mode_normalized = "stochastic_annual"
        elif mode_normalized == "intraday":
            mode_normalized = "stochastic_intraday"
            
    display_variables = normalize_display_variables(variables)
    
    # Load models
    # If stochastic mode, we need T, Q, P and Stochastic
    # If normal mode, we need T, Q, P (one of them is base_payload)
    
    if "stochastic" in mode_normalized:
        if not is_stochastic_file:
             # Try to find stochastic model if not passed directly
             # Assuming model_path is one of the harmonic models
             stem = model_path.stem
             # Remove suffix
             base_stem = stem
             for suffix in _TARGET_SUFFIXES.values():
                 if stem.endswith(suffix):
                     base_stem = stem[:-len(suffix)]
                     break
             stochastic_path = model_path.parent / f"{base_stem}_stochastic.json"
             if not stochastic_path.exists():
                 raise FileNotFoundError(f"Stochastic model not found at {stochastic_path} for mode {mode_normalized}")
        else:
             stochastic_path = model_path
             
        stochastic_payload = load_stochastic_model(stochastic_path)
        base_payload = stochastic_payload # Use this for metadata
    else:
        base_payload = load_linear_model(model_path)

    metadata = base_payload.get("metadata", {})
    station_name = str(
        metadata.get("station_usual_name")
        or metadata.get("station_name")
        or metadata.get("station_code")
        or ""
    ).strip()
    if not station_name:
        station_name = _slug_from_model_path(model_path).replace("_", " ").title()

    station_slug = slugify_station_name(station_name)
    if not station_slug:
        station_slug = _slug_from_model_path(model_path)

    artifact_paths = build_artifact_paths(station_slug)
    station_basename = compute_output_basename(station_slug)

    # Determine target variable for path override logic (less relevant for stochastic but good for consistency)
    target_variable_raw = str(metadata.get("target_variable") or "").strip().upper()
    if target_variable_raw in _TARGET_SUFFIXES:
        target_variable = target_variable_raw
    else:
        detected = None
        for candidate_target, candidate_suffix in _TARGET_SUFFIXES.items():
            if model_path.stem.endswith(candidate_suffix):
                detected = candidate_target
                break
        target_variable = detected or "T"

    if not is_stochastic_file:
        artifact_paths = _apply_model_path_override(artifact_paths, model_path, target_variable)

    model_map = {
        "T": artifact_paths.model_temperature_json,
        "Q": artifact_paths.model_specific_humidity_json,
        "P": artifact_paths.model_pressure_json,
    }

    payloads = {}
    # Load harmonic models
    for key, path in model_map.items():
        if not path.exists():
            raise FileNotFoundError(
                f"Missing required model JSON for target '{key}' at {path}. "
                "Run the generate pipeline first."
            )
        payloads[key] = load_linear_model(path)

    MEDIA_DIR.mkdir(parents=True, exist_ok=True)

    hist_daily = None
    hourly_hist = None
    real_df = None
    
    if artifact_paths.parquet.exists():
        try:
            sample_df = load_history_from_sample_data(artifact_paths.parquet)
            real_df = sample_df # Keep full df for stochastic comparison
            hist_daily, hourly_hist = historical_climatology_daily(sample_df)
        except Exception as exc:  # pragma: no cover - defensive logging
            print(f"[Warn] Failed to load historical dataset for overlays: {exc}")
            hist_daily = None
            hourly_hist = None
            real_df = None
    else:
        print(
            f"[Info] No cached dataset found at {artifact_paths.parquet}; plots will omit historical overlays.",
        )

    if mode_normalized == "intraday":
        media_path = MEDIA_DIR / f"{station_basename}_intraday_{day_int:03d}.png"
        plot_intraday(
            payloads["T"],
            payloads["Q"],
            payloads["P"],
            hourly_hist,
            station_name,
            day=day_int,
            save_path=media_path,
            variables=display_variables,
        )
    elif mode_normalized == "annual":
        days, stats = model_daily_stats_one_year_factorized(
            payloads["T"],
            payloads["Q"],
            payloads["P"],
        )

        media_path = MEDIA_DIR / f"{station_basename}_annual.png"
        plot_year(
            days,
            stats,
            station_name,
            media_path,
            hist_daily=hist_daily,
            variables=display_variables,
        )
    elif mode_normalized == "stochastic_annual":
        print(f"[Info] Generating stochastic annual simulation for {station_name}...")
        sim_df = simulate_stochastic_year(
            payloads["T"],
            payloads["Q"],
            payloads["P"],
            stochastic_payload,
            year=2010 # Use a complete historical year
        )
        media_path = MEDIA_DIR / f"{station_basename}_stochastic_annual.png"
        plot_stochastic_year(
            sim_df,
            station_name,
            save_path=media_path,
            real_df=real_df,
            variables=display_variables,
        )
    elif mode_normalized == "stochastic_intraday":
        print(f"[Info] Generating stochastic intraday simulation for {station_name} (Day {day_int})...")
        sim_df = simulate_stochastic_year(
            payloads["T"],
            payloads["Q"],
            payloads["P"],
            stochastic_payload,
            year=2010
        )
        media_path = MEDIA_DIR / f"{station_basename}_stochastic_intraday_{day_int:03d}.png"
        plot_stochastic_intraday(
            sim_df,
            station_name,
            day=day_int,
            save_path=media_path,
            real_df=real_df,
            variables=display_variables,
        )

    return media_path


def template_pipeline(model_name: str, target_language: str) -> Path:
    """Generate an embedded template for an existing model bundle."""

    ensure_directories()

    basename = _normalize_model_basename(model_name)
    language = target_language.strip().lower()

    if language != "cpp":
        raise ValueError(f"Unsupported target language '{target_language}'. Only 'cpp' is available.")

    temperature_path = MODEL_DIR / f"{basename}{_TARGET_SUFFIXES['T']}.json"
    specific_humidity_path = MODEL_DIR / f"{basename}{_TARGET_SUFFIXES['Q']}.json"
    pressure_path = MODEL_DIR / f"{basename}{_TARGET_SUFFIXES['P']}.json"

    if not temperature_path.exists():
        raise FileNotFoundError(f"Missing temperature model JSON at {temperature_path}.")
    if not specific_humidity_path.exists():
        raise FileNotFoundError(f"Missing specific humidity model JSON at {specific_humidity_path}.")
    if not pressure_path.exists():
        raise FileNotFoundError(f"Missing pressure model JSON at {pressure_path}.")

    temperature_payload = load_linear_model(temperature_path)
    specific_humidity_payload = load_linear_model(specific_humidity_path)
    pressure_payload = load_linear_model(pressure_path)

    stochastic_path = MODEL_DIR / f"{basename}_stochastic.json"
    stochastic_payload = None
    if stochastic_path.exists():
        try:
            stochastic_payload = load_stochastic_model(stochastic_path)
        except Exception as e:
            print(f"[Warn] Found stochastic model but failed to load it: {e}")
    else:
        print(f"[Info] No stochastic model found at {stochastic_path}; template will exclude it.")

    header_path = TEMPLATE_DIR / f"{basename}.hpp"
    generate_cpp_header(
        temperature_payload,
        specific_humidity_payload,
        pressure_payload,
        header_path,
        stochastic_payload=stochastic_payload,
    )

    return header_path


__all__ = [
    "ensure_directories",
    "clean_pipeline",
    "run_pipeline",
    "generate_pipeline",
    "regenerate_pipeline",
    "display_pipeline",
    "template_pipeline",
]
