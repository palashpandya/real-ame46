#!/usr/bin/env bash
# Run verification scripts 01-07 and 09-15 (08 is the slow stochastic search demo) and keep their output in results/.  Usage: ./run_all.sh [python]
set -u
PY=${1:-python3}
cd "$(dirname "$0")"
mkdir -p results
"$PY" - > results/environment.txt <<'EOF'
import platform, sys
print("# versions behind the logs in results/ (written by run_all.sh)")
print("python", platform.python_version(), platform.system().lower(), platform.machine())
for name, module in (("numpy", "numpy"), ("scipy", "scipy"), ("python-flint", "flint")):
    try:
        print(name, __import__(module).__version__)
    except ImportError:
        print(name, "not installed")
EOF
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-4}
status=0
for s in scripts/0[1-7]_*.py scripts/09_*.py scripts/1[0-5]_*.py; do
  name=$(basename "$s" .py)
  printf '%-28s' "$name"
  if "$PY" "$s" > "results/$name.log" 2>&1; then echo "ok"; else echo "FAILED (see results/$name.log)"; status=1; fi
done
exit $status
