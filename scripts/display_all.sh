#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODELS_DIR="${ROOT_DIR}/generated/models"

if command -v python >/dev/null 2>&1; then
  PYTHON_BIN="python"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="python3"
else
  echo "No Python interpreter found (expected 'python' or 'python3')." >&2
  exit 1
fi

if [[ ! -d "${MODELS_DIR}" ]]; then
  echo "No generated models directory found at ${MODELS_DIR}" >&2
  exit 1
fi

shopt -s nullglob
MODEL_FILES=("${MODELS_DIR}"/*.json)
shopt -u nullglob

if [[ ${#MODEL_FILES[@]} -eq 0 ]]; then
  echo "No model JSON files found in ${MODELS_DIR}" >&2
  exit 0
fi

declare -A STATION_ANNUAL=()
declare -A STATION_TEMPERATURE=()
declare -A STATION_STOCHASTIC=()

for model_path in "${MODEL_FILES[@]}"; do
  model_name="$(basename "${model_path}")"
  
  if [[ "${model_name}" == *"_stochastic.json" ]]; then
    station_key="${model_name%_stochastic.json}"
    STATION_STOCHASTIC["${station_key}"]="${model_path}"
    continue
  fi

  station_key="${model_name%.json}"

  if [[ "${model_name}" == *_temperature.json ]]; then
    station_key="${model_name%_temperature.json}"
    STATION_TEMPERATURE["${station_key}"]="${model_path}"
    STATION_ANNUAL["${station_key}"]="${model_path}"
    continue
  elif [[ "${model_name}" == *_specific_humidity.json ]]; then
    station_key="${model_name%_specific_humidity.json}"
  elif [[ "${model_name}" == *_pressure.json ]]; then
    station_key="${model_name%_pressure.json}"
  fi

  if [[ -z "${STATION_ANNUAL[${station_key}]+set}" ]]; then
    STATION_ANNUAL["${station_key}"]="${model_path}"
  fi
done

if [[ ${#STATION_ANNUAL[@]} -eq 0 ]]; then
  echo "No model JSON files found in ${MODELS_DIR}" >&2
  exit 0
fi

while IFS= read -r station; do
  model_path="${STATION_ANNUAL[${station}]}"
  model_name="$(basename "${model_path}")"
  echo "[HarmoClimate] Rendering annual plot for ${station} (${model_name})"
  "${PYTHON_BIN}" "${ROOT_DIR}/main.py" display "${model_path}"

  temp_path="${STATION_TEMPERATURE[${station}]:-}"
  if [[ -n "${temp_path}" ]]; then
    temp_name="$(basename "${temp_path}")"
    echo "[HarmoClimate] Rendering intraday plot (day 100) for ${station} (${temp_name})"
    "${PYTHON_BIN}" "${ROOT_DIR}/main.py" display "${temp_path}" --mode intraday --day 100
  else
    echo "[HarmoClimate] Skipping intraday plot for ${station} (no temperature bundle found)"
  fi

  stochastic_path="${STATION_STOCHASTIC[${station}]:-}"
  if [[ -n "${stochastic_path}" ]]; then
    stochastic_name="$(basename "${stochastic_path}")"
    echo "[HarmoClimate] Rendering stochastic annual plot for ${station} (${stochastic_name})"
    "${PYTHON_BIN}" "${ROOT_DIR}/main.py" display "${stochastic_path}" --mode stochastic_annual

    echo "[HarmoClimate] Rendering stochastic intraday plot (day 100) for ${station} (${stochastic_name})"
    "${PYTHON_BIN}" "${ROOT_DIR}/main.py" display "${stochastic_path}" --mode stochastic_intraday --day 100
  fi
done < <(printf "%s\n" "${!STATION_ANNUAL[@]}" | sort)
