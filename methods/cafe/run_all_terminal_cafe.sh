#!/usr/bin/env bash
set -euo pipefail

# Code follows this file; ROOT continues to identify the external data project.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

ROOT=${1:?Provide the external CAFE data project with inputs/cafe_counts and inputs/trees}
THREADS_PER_JOB=${THREADS_PER_JOB:-30}
PYTHON=${PYTHON:-python3}
CAFE=${CAFE:-cafe5}
OUT=${RESULT_ROOT:-"$ROOT/results"}
LOG=${LOG_ROOT:-"$ROOT/logs"}
RUN_FOCAL_AUDIT=${RUN_FOCAL_AUDIT:-1}

mkdir -p "$OUT" "$LOG"
exec > >(tee -a "$LOG/run_all_terminal_cafe.log") 2>&1
echo "run_start=$(date --iso-8601=seconds)"
echo "root=$ROOT"
echo "threads_per_job=$THREADS_PER_JOB"
echo "cafe=$($CAFE --version 2>&1 || true)"
echo "python=$($PYTHON --version 2>&1)"

run_model() {
  local dataset=$1
  local model=$2
  local tree=$3
  shift 3
  local model_out="$OUT/$dataset/$model"
  mkdir -p "$model_out"
  echo "[$dataset/$model] start=$(date --iso-8601=seconds)"
  "$CAFE" \
    -i "$ROOT/inputs/cafe_counts/${dataset}_counts.tsv" \
    -t "$ROOT/inputs/trees/${tree}.nwk" \
    -c "$THREADS_PER_JOB" -o "$model_out/results" -P 0.05 -p "$@"
  "$PYTHON" "$SCRIPT_DIR/parse_cafe_all_branches.py" \
    --model-name "$model" --results-dir "$model_out/results" --out-dir "$model_out/parsed"
  echo "[$dataset/$model] finish=$(date --iso-8601=seconds)"
}

run_dataset() {
  local dataset=$1
  local dataset_log="$LOG/${dataset}.log"
  {
    run_model "$dataset" base_corr_l1 correlated_l1
    run_model "$dataset" error_corr_l1 correlated_l1 -e
    local error_file
    error_file=$(find "$OUT/$dataset/error_corr_l1/results" -maxdepth 1 -type f -name '*_error_model.txt' | head -n 1)
    [[ -n "$error_file" ]] || { echo "Missing error model for $dataset" >&2; exit 3; }
    run_model "$dataset" gamma3_corr_l1 correlated_l1 -k 3
    run_model "$dataset" gamma3_error_corr_l1 correlated_l1 -k 3 "-e${error_file}"
    run_model "$dataset" base_corr_l01 correlated_l01
    run_model "$dataset" base_relaxed_l1 relaxed_l1
    "$PYTHON" "$SCRIPT_DIR/compare_terminal_models.py" \
      --dataset-dir "$OUT/$dataset" --out-dir "$OUT/$dataset/summary"
  } 2>&1 | tee -a "$dataset_log"
}

run_dataset original &
pid_original=$!
run_dataset augmented &
pid_augmented=$!
wait "$pid_original"
wait "$pid_augmented"

if [[ "$RUN_FOCAL_AUDIT" == "1" ]]; then
  "$PYTHON" "$SCRIPT_DIR/verify_focal_regression.py" --root "$ROOT" --results-root "$OUT"
fi
touch "$OUT/all_terminal_cafe.complete"
echo "run_finish=$(date --iso-8601=seconds)"
df -h "$ROOT" || true
