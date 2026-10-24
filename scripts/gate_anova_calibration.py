#!/usr/bin/env python3
"""S4.1 gates R-D (calibration) and R-F (power) of docs/rank.md section 4 for `stats test anova`:

    python3 scripts/gate_anova_calibration.py [--bin build/stats] [--replicates 4000] [--seed 1] [--workers 8]

Under a true null (k = 3 and 5 normal groups of mean 0, seeded), for sizes and standard deviations that differ between the groups, the rate at which `p < 0.05` lies inside the two-sided 99.9% binomial interval for the
number of replicates and the p-values are uniform (Kolmogorov-Smirnov, p > 0.001): the classic F where the variances are equal (it is an exact test there), Welch's F* on every configuration with at least 10 values in every group
(its approximation is printed, not gated, below). CONTROLS from the same replicates must be refused by the same test at some configuration, so the gate is able to fail: the classic F where the variances are not equal (and
the sizes are not), the degrees of freedom swapped, d2 + 1, the chi-square approximation of F d1. Power: equal variances and planted means, the rate at which `p < 0.05` of the classic F equals the power of the noncentral F (scipy,
an independent reference), inside the binomial interval.

Exit 0: every test passes and every control is refused. Exit 1: a test fails. Exit 2: the harness failed.
"""
import argparse
import concurrent.futures
import json
import pathlib
import subprocess
import sys
import tempfile

import scipy.stats as ss

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
from gate_test_calibration import ALPHA_TEST, interval, normal  # noqa: E402

# (sizes, sds): the classic F is judged where the sds are equal, Welch's on every configuration with 10 or more values in every group
CONFIGS = [
    ((10, 10, 10), (1.0, 1.0, 1.0)),
    ((5, 15, 40), (1.0, 1.0, 1.0)),
    ((4, 8, 16, 32, 64), (1.0, 1.0, 1.0, 1.0, 1.0)),
    ((10, 20, 40), (3.0, 2.0, 1.0)),
    ((40, 20, 10), (3.0, 2.0, 1.0)),
    ((12, 12, 12, 12, 12), (1.0, 2.0, 3.0, 4.0, 5.0)),
    ((10, 30, 50, 70, 90), (5.0, 4.0, 3.0, 2.0, 1.0)),
    ((2, 2, 2), (1.0, 1.0, 1.0)),
    ((3, 3, 3), (1.0, 1.0, 1.0)),
    ((5, 5, 5), (1.0, 3.0, 1.0)),
]


def write(d, groups):
    rows = ["x,g"]
    for gi, vals in enumerate(groups):
        rows += ["%r,g%d" % (v, gi) for v in vals]
    (d / "f.csv").write_text("\n".join(rows) + "\n")


def call(binary, d, variance):
    p = subprocess.run([binary, "test", "anova", "--root", str(d), "--cols", "x", "--group", "g", "--variance", variance, "f.csv"], capture_output=True, text=True, timeout=120)
    if p.returncode != 0:
        raise RuntimeError("exit %d: %s" % (p.returncode, p.stdout[:200]))
    return json.loads(p.stdout)["data"]["tests"][0]


def one(binary, d, groups):
    write(d, groups)
    return call(binary, d, "welch"), call(binary, d, "pooled")


def run_many(binary, jobs, workers):
    dirs = [pathlib.Path(tempfile.mkdtemp()) for _ in range(workers)]
    try:
        def work(k):
            return [one(binary, dirs[k], jobs[i]) for i in range(k, len(jobs), workers)]
        with concurrent.futures.ThreadPoolExecutor(workers) as ex:
            parts = list(ex.map(work, range(workers)))
        out = [None] * len(jobs)
        for k, part in enumerate(parts):
            for j, v in enumerate(part):
                out[k + j * workers] = v
        return out
    finally:
        for d in dirs:
            for f in d.glob("*"):
                f.unlink()
            d.rmdir()


def judge(ps, lo, hi):
    rej = sum(1 for p in ps if p < 0.05)
    ks = ss.kstest(ps, "uniform").pvalue
    return lo <= rej <= hi and ks > ALPHA_TEST, rej, ks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--replicates", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    reps = a.replicates
    lo, hi = interval(reps, 0.05)
    bad = []
    lines = []
    refused = {"classic-for-unequal-variances": False, "df-swapped": False, "df2-plus-one": False, "chi-square-approximation": False}
    try:
        for ci, (sizes, sds) in enumerate(CONFIGS):
            r = SM64(a.seed * 1000 + ci)
            jobs = [[[s * normal(r) for _ in range(n)] for n, s in zip(sizes, sds)] for _ in range(reps)]
            res = run_many(a.bin, jobs, a.workers)
            k = len(sizes)
            equal = len(set(sds)) == 1
            wp = [x[0]["p_value"] for x in res]
            cp = [x[1]["p_value"] for x in res]
            okw, rejw, ksw = judge(wp, lo, hi)
            gated = min(sizes) >= 10
            lines.append("Welch n=%s sd=%s: rejected %d of %d (interval %d..%d), KS p %.3f%s" % (sizes, sds, rejw, reps, lo, hi, ksw, ("" if okw else "  <-- FAILS") if gated else ("" if okw else "  (printed, not gated: a group below 10)")))
            if not okw and gated:
                bad.append("Welch n=%s sd=%s null calibration" % (sizes, sds))
            okc, rejc, ksc = judge(cp, lo, hi)
            if equal:
                lines.append("classic n=%s sd=%s: rejected %d of %d, KS p %.3f%s" % (sizes, sds, rejc, reps, ksc, "" if okc else "  <-- FAILS"))
                if not okc:
                    bad.append("classic n=%s null calibration" % (sizes,))
            elif not okc:
                refused["classic-for-unequal-variances"] = True
                lines.append("classic n=%s sd=%s: rejected %d of %d, KS p %.3f (outside, as it must be: its assumption fails)" % (sizes, sds, rejc, reps, ksc))
            # controls: the classic replicates (an exact F) read wrongly
            fs = [x[1]["statistic"] for x in res]
            d2 = sum(sizes) - k
            if equal:
                ok_s, _, _ = judge([float(ss.f.sf(f, d2, k - 1)) for f in fs], lo, hi)
                ok_p, _, _ = judge([float(ss.f.sf(f, k - 1, d2 + 1)) for f in fs], lo, hi)
                ok_c, _, _ = judge([float(ss.chi2.sf(f * (k - 1), k - 1)) for f in fs], lo, hi)
                if not ok_s and k - 1 != d2:
                    refused["df-swapped"] = True
                if not ok_p:
                    refused["df2-plus-one"] = True
                if not ok_c:
                    refused["chi-square-approximation"] = True
        # power: equal variances, planted means, the classic F against the noncentral F
        sizes, means = (15, 20, 25), (0.0, 0.5, 1.0)
        pr = SM64(a.seed * 7919)
        jobs = [[[m + normal(pr) for _ in range(n)] for n, m in zip(sizes, means)] for _ in range(reps)]
        res = run_many(a.bin, jobs, a.workers)
        n_tot = sum(sizes)
        grand = sum(n * m for n, m in zip(sizes, means)) / n_tot
        lam = sum(n * (m - grand) ** 2 for n, m in zip(sizes, means))
        d1, d2 = len(sizes) - 1, n_tot - len(sizes)
        power = float(ss.ncf.sf(ss.f.isf(0.05, d1, d2), d1, d2, lam))
        rej = sum(1 for x in res if x[1]["p_value"] < 0.05)
        plo, phi = interval(reps, power)
        ok = plo <= rej <= phi
        lines.append("planted means (0, 0.5, 1), n=(15,20,25), classic: power (noncentral F, scipy) %.6f, rejected %d of %d (interval %d..%d)%s" % (power, rej, reps, plo, phi, "" if ok else "  <-- FAILS"))
        if not ok:
            bad.append("power of the planted means")
    except Exception as e:
        print("gate_anova_calibration: harness failure: %s" % e)
        return 2
    for name, was in refused.items():
        lines.append("control %s: %s" % (name, "refused by the calibration test" if was else "NOT refused"))
        if not was:
            bad.append("control %s was not refused: the gate cannot tell it from the program" % name)
    print("gate_anova_calibration (%d replicates per configuration, seed %d):" % (reps, a.seed))
    for ln in lines:
        print("  " + ln)
    print("gate_anova_calibration: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
