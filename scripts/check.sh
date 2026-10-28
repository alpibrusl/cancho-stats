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
step "stats test"; $PY scripts/gate_test.py --files 200 --scipy; $PY scripts/gate_test_calibration.py
step "stats test t2 and paired"; $PY scripts/gate_test2.py --kind t2 --files 200 --scipy; $PY scripts/gate_test2.py --kind paired --files 200 --scipy; $PY scripts/gate_test2_calibration.py --kind t2; $PY scripts/gate_test2_calibration.py --kind paired
step "stats test chisq"; $PY scripts/gate_test3.py --files 300 --scipy; $PY scripts/gate_test3_calibration.py
step "stats test prop"; $PY scripts/gate_test4.py --files 300; $PY scripts/gate_test4_calibration.py
step "stats test --correct"; $PY scripts/gate_test5.py --files 100; $PY scripts/gate_test5_calibration.py --replicates 1000
step "stats test anova"; $PY scripts/gate_anova.py --files 200 --scipy; $PY scripts/gate_anova_calibration.py --replicates 1000
step "stats test kruskal"; $PY scripts/gate_kruskal.py --files 200 --scipy; $PY scripts/gate_kruskal_calibration.py --replicates 1000
step "stats ci prop"; $PY scripts/gate_ci_prop.py --files 200; $PY scripts/gate_ci_prop_coverage.py; $PY scripts/gate_ci_mean.py --files 200
step "stats corr --exact"; $PY scripts/gate_corr_exact.py --files 200 --scipy
step "stats test fisher"; $PY scripts/gate_fisher.py --files 200 --scipy
step "stats test ks"; $PY scripts/gate_ks.py --files 200 --scipy; $PY scripts/gate_ks_calibration.py --replicates 1000
step "stats test mannwhitney"; $PY scripts/gate_mannwhitney.py --files 200 --scipy; $PY scripts/gate_mannwhitney_calibration.py --replicates 1000
step "stats test wilcoxon"; $PY scripts/gate_wilcoxon.py --files 200 --scipy; $PY scripts/gate_wilcoxon_calibration.py --replicates 1000
step "mutant patterns"; for m in bigint psum xsum stats central quantile histogram levels by cols corr test test2 test3 test4 test5 anova kruskal mannwhitney wilcoxon ks fisher corrx ci special; do $PY scripts/mutants_$m.py --check; done
step "exact accumulation"; $PY scripts/measure_acc.py | tail -3
step "accuracy sweeps (fast)"; $PY scripts/gate_special.py --fast
step "independence"; $PY scripts/gate_independence.py --fast
step "determinism"; $PY scripts/xplat_fp.py; $PY scripts/xplat_fp.py --mutants | tail -2
if [ "$1" = "--full" ]; then
  step "accuracy sweeps (10%)"; $PY scripts/gate_special.py --n-scale 0.1
  step "mutants of the special functions"; $PY scripts/mutants_special.py
  step "mutants of the big integers, the power sums and stats"; $PY scripts/mutants_bigint.py; $PY scripts/mutants_psum.py; $PY scripts/mutants_xsum.py; $PY scripts/mutants_stats.py; $PY scripts/mutants_central.py; $PY scripts/mutants_quantile.py; $PY scripts/mutants_histogram.py; $PY scripts/mutants_levels.py; $PY scripts/mutants_by.py; $PY scripts/mutants_cols.py; $PY scripts/mutants_corr.py; $PY scripts/mutants_test.py; $PY scripts/mutants_test2.py; $PY scripts/mutants_test3.py; $PY scripts/mutants_test4.py; $PY scripts/mutants_test5.py; $PY scripts/mutants_anova.py; $PY scripts/mutants_kruskal.py; $PY scripts/mutants_mannwhitney.py; $PY scripts/mutants_wilcoxon.py; $PY scripts/mutants_ks.py; $PY scripts/mutants_fisher.py; $PY scripts/mutants_corrx.py; $PY scripts/mutants_ci.py
fi
echo; echo "check.sh: all gates passed"
