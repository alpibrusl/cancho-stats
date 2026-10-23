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
step "build"; $CANCHO install >/dev/null; $CANCHO build --bin special-check; $CANCHO build --bin acc-check; $CANCHO build --bin bigint-check; $CANCHO build --bin psum-check; $CANCHO build --bin xsum-check; $CANCHO build --bin stats; $CANCHO build --bin quantile-check
step "consumption"; $PY scripts/consume.py --mutants | tail -2
step "domain and fuzz"; $PY scripts/gate_domain.py --fuzz 50000 | grep -v "first case with"
step "quantile properties"; $PY scripts/gate_qprops.py --n 2000
step "big integers"; $PY scripts/gate_bigint.py --n 3000
step "quantiles"; $PY scripts/gate_quantile.py --cases 3000
step "stats describe"; $PY scripts/gate_describe.py --files 200; $PY scripts/gate_describe_refusals.py --fuzz 1000; $PY scripts/manifest.py --check; $PY scripts/update_replay.py --check
step "power sums"; $PY scripts/gate_psum.py --vectors 300; $PY scripts/gate_xsum.py --vectors 300
step "stats corr"; $PY scripts/gate_corr.py --files 200 --scipy; $PY scripts/gate_corr.py --files 200 --scipy --method spearman; $PY scripts/gate_corr.py --files 200 --scipy --method kendall; $PY scripts/gate_corr_calibration.py; $PY scripts/gate_corr_calibration.py --method spearman; $PY scripts/gate_corr_calibration.py --method kendall
step "exact accumulation"; $PY scripts/measure_acc.py | tail -3
step "accuracy sweeps (fast)"; $PY scripts/gate_special.py --fast
step "independence"; $PY scripts/gate_independence.py --fast
step "determinism"; $PY scripts/xplat_fp.py; $PY scripts/xplat_fp.py --mutants | tail -2
if [ "$1" = "--full" ]; then
  step "accuracy sweeps (10%)"; $PY scripts/gate_special.py --n-scale 0.1
  step "mutants of the special functions"; $PY scripts/mutants_special.py
  step "mutants of the big integers, the power sums and stats"; $PY scripts/mutants_bigint.py; $PY scripts/mutants_psum.py; $PY scripts/mutants_xsum.py; $PY scripts/mutants_stats.py; $PY scripts/mutants_central.py; $PY scripts/mutants_quantile.py; $PY scripts/mutants_histogram.py; $PY scripts/mutants_levels.py; $PY scripts/mutants_by.py; $PY scripts/mutants_cols.py; $PY scripts/mutants_corr.py
fi
echo; echo "check.sh: all gates passed"
