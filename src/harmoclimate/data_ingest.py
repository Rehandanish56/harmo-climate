"""Data ingestion utilities for HarmoClimate."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError as exc:  # pragma: no cover - execution stops before tests
    raise ImportError("pyarrow is required for HarmoClimate data ingestion.") from exc

from .config import CHUNK_SIZE, STATION_CODE, build_artifact_paths, slugify_station_name
from .core import DATASET_COLUMNS


@dataclass
class StationRecord:
    """Lightweight station-level record captured while streaming."""

    station_code: str | None
    station_name: str
    lon: float
    lat: float
    alti: float
    delta_utc_solar_h: float


@dataclass
class StreamResult:
    """Container bundling streamed station records and generated artefact info."""

    station_records: list[StationRecord]
    parquet_path: Path
    station_name: str
    station_slug: str


def choose_columns(colnames: Iterable[str]) -> tuple[str, str, str, str, str, str, str, str, str]:
    """Validate the expected Meteo-France column names."""

    col_code = "NUM_POSTE"
    col_station = "NOM_USUEL"
    col_dt = "AAAAMMJJHH"
    col_T = "T"
    col_U = "U"
    col_lon = "LON"
    col_lat = "LAT"
    col_alti = "ALTI"
    col_p = "PSTAT"

    required = [col_code, col_station, col_dt, col_T, col_U, col_lon, col_lat, col_alti, col_p]
    missing = [c for c in required if c not in colnames]
    if missing:
        raise RuntimeError(f"Missing required columns: {missing} in {colnames}")

    return col_code, col_station, col_dt, col_T, col_U, col_lon, col_lat, col_alti, col_p


def parse_dt_aaaammjjhh(series: pd.Series) -> pd.Series:
    """Parse AAAAMMJJHH timestamps expressed in local French time."""

    # Coerce to string and handle potential float representation (e.g. 2021010100.0)
    s = series.astype(str).str.strip()
    
    # Remove trailing .0 if present (common in CSVs loaded as float)
    s = s.str.replace(r"\.0$", "", regex=True)
    
    # AAAAMMJJHH is 10 chars.
    s = s.str.zfill(10)

    return pd.to_datetime(s, format="%Y%m%d%H", errors="coerce")


def _normalize_station_code(value: object) -> str:
    """Coerce a station code value to the canonical string representation."""

    code = str(value).strip()
    if not code:
        return ""
    code = code.replace(".0", "")  # CSV exports frequently store codes as floats
    if code.lower() in {"nan", "na", "<na>", "none"}:
        return ""
    return code


def _process_chunk(
    chunk: pd.DataFrame,
    station_code_str: str,
) -> tuple[pd.DataFrame | None, list[StationRecord]]:
    """
    Process a single chunk of data: filter by station, clean types, and extract records.
    
    Returns:
        Tuple of (processed_dataframe, list_of_station_records).
        processed_dataframe is None if no data remains after filtering.
    """
    (
        col_code,
        col_station,
        col_dt,
        col_T,
        col_U,
        col_lon,
        col_lat,
        col_alti,
        col_p,
    ) = choose_columns(list(chunk.columns))

    # Filter by station code

    code_series = chunk[col_code].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
    mask = code_series == station_code_str
    
    if not mask.any():
        return None, []

    keep_cols = [col_dt, col_T, col_U, col_lon, col_station, col_lat, col_code, col_alti, col_p]
    sub = chunk.loc[mask, keep_cols].copy()
    
    if sub.empty:
        return None, []

    # Type conversion and cleaning
    sub["dt_local"] = parse_dt_aaaammjjhh(sub[col_dt])
    

    numeric_cols = {
        "T": col_T,
        "RH": col_U,
        "LON": col_lon,
        "LAT": col_lat,
        "ALTI": col_alti,
        "P": col_p
    }
    
    for target, source in numeric_cols.items():
        sub[target] = pd.to_numeric(sub[source], errors="coerce")

    sub["STATION_CODE"] = code_series[mask].values
    sub["STATION_NAME"] = sub[col_station].astype(str).str.strip()

    # Drop invalid rows
    required_cols = ["dt_local", "T", "RH", "LON", "LAT", "ALTI", "P"]
    sub = sub.dropna(subset=required_cols)
    
    if sub.empty:
        return None, []

    # Timezone handling
    sub["dt_local"] = sub["dt_local"].dt.tz_localize(
        "Europe/Paris",
        nonexistent="shift_forward",
        ambiguous="NaT",
    )
    sub = sub.dropna(subset=["dt_local"])
    
    if sub.empty:
        return None, []

    sub["DT_UTC"] = sub["dt_local"].dt.tz_convert("UTC")

    # Final type casting
    sub["T"] = sub["T"].astype("float32")
    sub["RH"] = sub["RH"].clip(0, 100).astype("float32")
    sub["P"] = sub["P"].astype("float32")
    sub["LON"] = sub["LON"].astype("float32")
    sub["LAT"] = sub["LAT"].astype("float32")
    sub["ALTI"] = sub["ALTI"].astype("float32")
    sub["STATION_CODE"] = sub["STATION_CODE"].astype("string")
    sub["STATION_NAME"] = sub["STATION_NAME"].astype("string")

    # Extract StationRecords
    
    records_df = pd.DataFrame({
        "station_code": sub["STATION_CODE"].astype(object), # keep as object for dataclass
        "station_name": sub["STATION_NAME"].astype(object),
        "lon": sub["LON"].astype(float),
        "lat": sub["LAT"].astype(float),
        "alti": sub["ALTI"].astype(float),
        "delta_utc_solar_h": (sub["LON"] / 15.0).astype(float)
    })
    
    # Handle potential NaNs in names for the record
    records_df["station_name"] = records_df["station_name"].fillna("")
    
    # Convert to list of dataclasses
    station_records = [StationRecord(**row) for row in records_df.to_dict("records")]

    # Select final columns for Parquet
    sub = sub[list(DATASET_COLUMNS)]
    
    return sub, station_records


def stream_filter_to_disk(
    urls: Iterable[str],
    station_code: str = STATION_CODE,
    chunk_size: int = CHUNK_SIZE,
) -> StreamResult:
    """Stream Meteo-France archives, filter rows, and persist the reduced dataset."""

    station_code_str = _normalize_station_code(station_code)

    writer: pq.ParquetWriter | None = None
    parquet_path: Path | None = None
    station_name: str | None = None
    station_slug: str | None = None

    total_kept = 0
    all_station_records: list[StationRecord] = []

    for url in urls:
        print(f"[Stream] {url}")
        reader = pd.read_csv(
            url,
            compression="gzip",
            sep=";",
            encoding="utf-8",
            chunksize=chunk_size,
            low_memory=True,
        )
        for i, chunk in enumerate(reader, 1):
            processed_df, records = _process_chunk(chunk, station_code_str)
            
            if processed_df is None or processed_df.empty:
                continue

            all_station_records.extend(records)

            # Initialize output info on first valid data found
            if station_name is None:
                # Try to find a valid name in the records
                candidate = next((r.station_name for r in records if r.station_name), "")
                fallback_name = station_code_str if not candidate else candidate
                station_name = candidate or fallback_name
                station_slug = slugify_station_name(station_name) or slugify_station_name(fallback_name)
                
                parquet_path = build_artifact_paths(station_slug).parquet
                parquet_path.parent.mkdir(parents=True, exist_ok=True)
                parquet_path.unlink(missing_ok=True)

            table = pa.Table.from_pandas(processed_df, preserve_index=False)
            
            if parquet_path is None or station_slug is None:
                 # Should be unreachable given the init block above
                raise RuntimeError("Parquet output path was not initialised.")
                
            if writer is None:
                writer = pq.ParquetWriter(str(parquet_path), table.schema)
            
            writer.write_table(table)

            total_kept += len(processed_df)
            if i % 50 == 0:
                print(f"  chunks: {i:4d} | kept rows: {total_kept:,}")

    if writer is not None:
        writer.close()

    if parquet_path is None or station_name is None or station_slug is None:
        raise RuntimeError(f"No data found for station code {station_code_str}")

    print(f"[OK] Wrote {total_kept:,} filtered rows -> {parquet_path}")
    return StreamResult(
        station_records=all_station_records,
        parquet_path=parquet_path,
        station_name=station_name,
        station_slug=station_slug,
    )


__all__ = [
    "StationRecord",
    "StreamResult",
    "choose_columns",
    "parse_dt_aaaammjjhh",
    "stream_filter_to_disk",
]
