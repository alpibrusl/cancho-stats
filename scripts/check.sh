#!/bin/sh
# Every gate of milestone S0 (docs/numerics.md section 7), in the order that fails fastest first. Needs the compiler ($CANCHO, default ../cancho-compiler-a4572ea/target/release/cancho),
# the pinned stores (`cancho install`) and the venv of ci/requirements.txt in .venv (python3 -m venv .venv && .venv/bin/pip install -r ci/requirements.txt).
#   sh scripts/check.sh          # the fast gates (minutes)
#   sh scripts/check.sh --full   # also the sweeps at 10% of the full size and the mutants of the special functions (an hour or more)
set -e
cd "$(dirname "$0")/.."
PY=${PY:-.venv/bin/python}
CANCHO=${CANCHO:-../cancho-compiler-a4572ea/target/release/cancho}
export CANCHO
step() { echo; echo "== $*"; }
step "pins"; sh scripts/check_pins.sh; $PY scripts/mutants_pins.py | tail -1
step "workflow"; $PY scripts/check_workflow.py; $PY scripts/check_workflow.py --mutants | tail -1
step "documents"; sh scripts/check_docs.sh
step "references"; $PY scripts/special_ref.py --selfcheck
step "build"; $CANCHO install >/dev/null; $CANCHO build --bin special-check; $CANCHO build --bin acc-check; $CANCHO build --bin bigint-check
step "consumption"; $PY scripts/consume.py --mutants | tail -2
step "domain and fuzz"; $PY scripts/gate_domain.py --fuzz 50000 | grep -v "first case with"
step "big integers"; $PY scripts/gate_bigint.py --n 3000
step "exact accumulation"; $PY scripts/measure_acc.py | tail -3
step "accuracy sweeps (fast)"; $PY scripts/gate_special.py --fast
step "independence"; $PY scripts/gate_independence.py --fast
step "determinism"; $PY scripts/xplat_fp.py; $PY scripts/xplat_fp.py --mutants | tail -2
if [ "$1" = "--full" ]; then
  step "accuracy sweeps (10%)"; $PY scripts/gate_special.py --n-scale 0.1
  step "mutants of the special functions"; $PY scripts/mutants_special.py
fi
echo; echo "check.sh: all gates passed"
