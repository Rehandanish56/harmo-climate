#!/usr/bin/env bash

set -euo pipefail

if command -v python >/dev/null 2>&1; then
  PYTHON_BIN="python"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="python3"
else
  echo "No Python interpreter found (expected 'python' or 'python3')." >&2
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODELS_DIR="${ROOT_DIR}/generated/models"

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

declare -A processed_stations=()

for model_path in "${MODEL_FILES[@]}"; do
  model_name="$(basename "${model_path}")"
  
  if [[ "${model_name}" == *"_stochastic.json" ]]; then
    continue
  fi

  stem="${model_name%.json}"

  base_slug="${stem}"
  for suffix in "_temperature" "_specific_humidity" "_pressure"; do
    if [[ "${base_slug}" == *"${suffix}" ]]; then
      base_slug="${base_slug%${suffix}}"
      break
    fi
  done

  if [[ -n "${processed_stations["${base_slug}"]+yes}" ]]; then
    continue
  fi

  processed_stations["${base_slug}"]=1

  echo "[HarmoClimate] Regenerating ${base_slug} (via ${model_name})"
  "${PYTHON_BIN}" "${ROOT_DIR}/main.py" regenerate "${model_path}"
done
