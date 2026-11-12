#!/usr/bin/env bash
set -euo pipefail

# Code follows this file; ROOT continues to identify the external data project.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

ROOT=${1:?Provide the external CAFE data project with inputs/cafe_counts and inputs/trees}
REPLICATE=${2:?Usage: run_multistart_replicate.sh ROOT replicate_XX}
ENV_PREFIX=${ENV_PREFIX:-$HOME/envs/fachb416_cafe5}
OUT="$ROOT/multistart/$REPLICATE/results"
LOG="$ROOT/multistart/$REPLICATE/logs"

if [[ -e "$OUT" ]] || [[ -e "$LOG" ]]; then
  echo "Refusing to overwrite existing replicate: $REPLICATE" >&2
  exit 3
fi
mkdir -p "$LOG"

export PYTHON="$ENV_PREFIX/bin/python"
export CAFE="$ENV_PREFIX/bin/cafe5"
export THREADS_PER_JOB=${THREADS_PER_JOB:-30}
export RESULT_ROOT="$OUT"
export LOG_ROOT="$LOG"
export RUN_FOCAL_AUDIT=0

echo "replicate_start=$(date --iso-8601=seconds)" | tee "$LOG/replicate_driver.log"
bash "$SCRIPT_DIR/run_all_terminal_cafe.sh" "$ROOT"
echo "replicate_finish=$(date --iso-8601=seconds)" | tee -a "$LOG/replicate_driver.log"
