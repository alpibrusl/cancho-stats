#!/usr/bin/env python3
"""S3.2 gates T-E (calibration) and T-F (power and coverage) of docs/test.md section 5 for `stats test t2` and `stats test paired`:

    python3 scripts/gate_test2_calibration.py [--bin build/stats] [--kind t2|paired] [--replicates 4000] [--seed 1] [--workers 8]

t2. Under a true null (two normal groups of mean 0, seeded), for sizes and standard deviations that differ between the groups, the rate at which `p < 0.05` of Welch's test lies inside the two-sided 99.9% binomial
interval for the number of replicates, the p-values are uniform (Kolmogorov-Smirnov, p > 0.001), and the 95% and 90% intervals cover the true difference at their rate (Welch's from 5 values per group on: its approximation is conservative below, printed). The pooled test is held to the same
where the variances are equal. CONTROLS from the same replicates must be refused by the same test at some configuration, so the gate is able to fail: the pooled p-value where Welch's is meant (unequal variances
and sizes), the degrees of freedom plus one, a one-sided p-value, the normal in place of Student's t. Power: with equal variances and a planted difference the rate at which `p < 0.05` equals the power of the
noncentral t (scipy, an independent reference), inside the binomial interval; the coverage of the planted difference is nominal for Welch's interval under unequal variances.
paired. Pairs that share a large common part (the usual reason to pair), n = 3, 5, 10, 30, 100, 1,000: the same tests as t1, and the control that the two columns were treated as independent samples must be refused.

Exit 0: every test passes and every control is refused. Exit 1: a test fails. Exit 2: the harness failed.
"""
import argparse
import concurrent.futures
import json
import math
import pathlib
import subprocess
import sys
import tempfile

import scipy.stats as ss

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
from gate_test_calibration import ALPHA_TEST, interval, normal  # noqa: E402

# (n1, n2, sd1, sd2): Welch is judged on every one, the pooled test where the sds are equal; unequal sizes with the larger variance in the smaller group, in the larger group, and in both
CONFIGS = [(2, 2, 1.0, 1.0), (3, 3, 1.0, 1.0), (5, 15, 1.0, 1.0), (30, 30, 1.0, 1.0), (10, 40, 4.0, 1.0), (40, 10, 4.0, 1.0), (20, 20, 1.0, 3.0), (100, 300, 2.0, 1.0)]
NS = (3, 5, 10, 30, 100, 1000)


def write_t2(d, a, b):
    rows = ["x,g"] + ["%r,a" % v for v in a] + ["%r,b" % v for v in b]
    (d / "f.csv").write_text("\n".join(rows) + "\n")


def call(binary, d, args):
    p = subprocess.run([binary, "test", "--root", str(d)] + args, capture_output=True, text=True, timeout=120)
    if p.returncode != 0:
        raise RuntimeError("exit %d: %s" % (p.returncode, p.stdout[:200]))
    return json.loads(p.stdout)["data"]["tests"][0]


def one_t2(binary, d, job):
    a, b = job
    write_t2(d, a, b)
    base = ["t2", "--cols", "x", "--group", "g"]
    return (call(binary, d, base + ["--variance", "welch", "--level", "0.95", "f.csv"]),
            call(binary, d, base + ["--variance", "welch", "--level", "0.9", "f.csv"]),
            call(binary, d, base + ["--variance", "pooled", "--level", "0.95", "f.csv"]))


def one_paired(binary, d, job):
    a, b = job
    (d / "f.csv").write_text("a,b\n" + "\n".join("%r,%r" % (x, y) for x, y in zip(a, b)) + "\n")
    return (call(binary, d, ["paired", "--pairs", "a:b", "--level", "0.95", "f.csv"]),
            call(binary, d, ["paired", "--pairs", "a:b", "--level", "0.9", "f.csv"]))


def run_many(fn, binary, jobs, workers):
    dirs = [pathlib.Path(tempfile.mkdtemp()) for _ in range(workers)]
    try:
        def work(k):
            return [fn(binary, dirs[k], jobs[i]) for i in range(k, len(jobs), workers)]
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


def calibrate_t2(a, lines, bad, refused):
    reps = a.replicates
    lo, hi = interval(reps, 0.05)
    clo, chi = interval(reps, 0.95)
    c9lo, c9hi = interval(reps, 0.9)
    for ci, (n1, n2, s1, s2) in enumerate(CONFIGS):
        r = SM64(a.seed * 1000 + ci)
        jobs = [([s1 * normal(r) for _ in range(n1)], [s2 * normal(r) for _ in range(n2)]) for _ in range(reps)]
        res = run_many(one_t2, a.bin, jobs, a.workers)
        w95 = [x[0] for x in res]
        ps = [x["p_value"] for x in w95]
        ok, rej, ks = judge(ps, lo, hi)
        cov95 = sum(1 for x in w95 if x["interval"]["low"] <= 0.0 <= x["interval"]["high"])
        cov90 = sum(1 for x in res if x[1]["interval"]["low"] <= 0.0 <= x[1]["interval"]["high"])
        ok = ok and clo <= cov95 <= chi and c9lo <= cov90 <= c9hi
        # Welch's approximation is itself conservative with a group of 3 (the rate of 3.2% at 5% in the first run, the p-value being the exact Welch formula, docs/test.md section 9): printed, not gated, below 5 per group
        gated = min(n1, n2) >= 5
        lines.append("Welch n=(%d,%d) sd=(%g,%g): rejected %d of %d (interval %d..%d), KS p %.3f; the 95%% interval covered %d (%d..%d), the 90%% interval %d (%d..%d)%s" % (n1, n2, s1, s2, rej, reps, lo, hi, ks, cov95, clo, chi, cov90, c9lo, c9hi, ("" if ok else "  <-- FAILS") if gated else ("" if ok else "  (printed, not gated: a group below 5)")))
        if not ok and gated:
            bad.append("Welch n=(%d,%d) sd=(%g,%g) null calibration or coverage" % (n1, n2, s1, s2))
        if s1 == s2:
            pp = [x[2]["p_value"] for x in res]
            okp, rejp, ksp = judge(pp, lo, hi)
            covp = sum(1 for x in res if x[2]["interval"]["low"] <= 0.0 <= x[2]["interval"]["high"])
            okp = okp and clo <= covp <= chi
            # the control of the degrees of freedom is judged where the p-value is an exact Student t: the pooled test, whose smallest case (3,3) has 4 (Welch's own p-value is only approximately t)
            okc, _, _ = judge([2.0 * ss.t.sf(abs(x[2]["statistic"]), x[2]["df"] + 1.0) for x in res], lo, hi)
            if not okc:
                refused["df-plus-one"] = True
            lines.append("pooled n=(%d,%d) sd=(%g,%g): rejected %d of %d, KS p %.3f; the 95%% interval covered %d%s" % (n1, n2, s1, s2, rejp, reps, ksp, covp, "" if okp else "  <-- FAILS"))
            if not okp:
                bad.append("pooled n=(%d,%d) null calibration or coverage" % (n1, n2))
        # controls: the same replicates read wrongly
        ts = [x["statistic"] for x in w95]
        dfs = [x["df"] for x in w95]
        wrong = {
            "pooled-for-welch": [x[2]["p_value"] for x in res],
            "df-plus-one": [2.0 * ss.t.sf(abs(t), d + 1.0) for t, d in zip(ts, dfs)],
            "one-sided": [p / 2.0 for p in ps],
            "normal": [2.0 * ss.norm.sf(abs(t)) for t in ts],
        }
        for name, cp in wrong.items():
            okc, _, _ = judge(cp, lo, hi)
            if not okc:
                refused[name] = True
    # power: equal variances, a planted difference, against the noncentral t
    n1, n2, delta = 20, 30, 0.6
    pr = SM64(a.seed * 7919)
    jobs = [([delta + normal(pr) for _ in range(n1)], [normal(pr) for _ in range(n2)]) for _ in range(reps)]
    res = run_many(one_t2, a.bin, jobs, a.workers)
    df = n1 + n2 - 2
    nc = delta / math.sqrt(1.0 / n1 + 1.0 / n2)
    tc = ss.t.ppf(0.975, df)
    power = float(ss.nct.sf(tc, df, nc) + ss.nct.cdf(-tc, df, nc))
    rej = sum(1 for x in res if x[2]["p_value"] < 0.05)
    plo, phi = interval(reps, power)
    cov = sum(1 for x in res if x[2]["interval"]["low"] <= delta <= x[2]["interval"]["high"])
    ok = plo <= rej <= phi and clo <= cov <= chi
    lines.append("planted difference 0.6 sd, n=(20,30), pooled: power (noncentral t, scipy) %.6f, rejected %d of %d (interval %d..%d); the 95%% interval covered the planted difference %d times (%d..%d)%s" % (power, rej, reps, plo, phi, cov, clo, chi, "" if ok else "  <-- FAILS"))
    if not ok:
        bad.append("power or coverage of the planted difference (pooled)")
    # coverage of a planted difference under unequal variances and sizes (Welch)
    n1, n2, s1, s2, delta = 10, 40, 4.0, 1.0, 2.0
    jobs = [([delta + s1 * normal(pr) for _ in range(n1)], [s2 * normal(pr) for _ in range(n2)]) for _ in range(reps)]
    res = run_many(one_t2, a.bin, jobs, a.workers)
    cov = sum(1 for x in res if x[0]["interval"]["low"] <= delta <= x[0]["interval"]["high"])
    ok = clo <= cov <= chi
    lines.append("planted difference 2, n=(10,40), sd=(4,1), Welch: the 95%% interval covered the planted difference %d times (%d..%d)%s" % (cov, clo, chi, "" if ok else "  <-- FAILS"))
    if not ok:
        bad.append("coverage of the planted difference (Welch, unequal variances)")


def calibrate_paired(a, lines, bad, refused):
    reps = a.replicates
    lo, hi = interval(reps, 0.05)
    clo, chi = interval(reps, 0.95)
    c9lo, c9hi = interval(reps, 0.9)
    for n in NS:
        r = SM64(a.seed * 1000 + n)
        jobs = []
        for _ in range(reps):
            z = [3.0 * normal(r) for _ in range(n)]
            jobs.append(([zi + normal(r) for zi in z], [zi + normal(r) for zi in z]))
        res = run_many(one_paired, a.bin, jobs, a.workers)
        ps = [x[0]["p_value"] for x in res]
        ts = [x[0]["statistic"] for x in res]
        ok, rej, ks = judge(ps, lo, hi)
        cov95 = sum(1 for x in res if x[0]["interval"]["low"] <= 0.0 <= x[0]["interval"]["high"])
        cov90 = sum(1 for x in res if x[1]["interval"]["low"] <= 0.0 <= x[1]["interval"]["high"])
        ok = ok and clo <= cov95 <= chi and c9lo <= cov90 <= c9hi
        lines.append("n=%d: rejected %d of %d (interval %d..%d), KS p %.3f; the 95%% interval covered %d (%d..%d), the 90%% interval %d (%d..%d)%s" % (n, rej, reps, lo, hi, ks, cov95, clo, chi, cov90, c9lo, c9hi, "" if ok else "  <-- FAILS"))
        if not ok:
            bad.append("n=%d null calibration or coverage" % n)
        # the two columns as independent samples: the difference is the same, the standard error is not
        indep = []
        for (xa, xb) in jobs:
            indep.append(float(ss.ttest_ind(xa, xb, equal_var=False).pvalue))
        wrong = {
            "independent-samples": indep,
            "df-plus-one": [2.0 * ss.t.sf(abs(t), n) for t in ts],
            "one-sided": [p / 2.0 for p in ps],
            "normal": [2.0 * ss.norm.sf(abs(t)) for t in ts],
        }
        for name, cp in wrong.items():
            okc, _, _ = judge(cp, lo, hi)
            if not okc:
                refused[name] = True
    n, delta = 30, 0.5
    pr = SM64(a.seed * 7919)
    jobs = []
    for _ in range(reps):
        z = [3.0 * normal(pr) for _ in range(n)]
        # the differences have sd sqrt(2): plant a mean difference of 0.5 sd of the differences
        d = delta * math.sqrt(2.0)
        jobs.append(([zi + normal(pr) + d for zi in z], [zi + normal(pr) for zi in z]))
    res = run_many(one_paired, a.bin, jobs, a.workers)
    tc = ss.t.ppf(0.975, n - 1)
    nc = delta * math.sqrt(n)
    power = float(ss.nct.sf(tc, n - 1, nc) + ss.nct.cdf(-tc, n - 1, nc))
    rej = sum(1 for x in res if x[0]["p_value"] < 0.05)
    plo, phi = interval(reps, power)
    d = delta * math.sqrt(2.0)
    cov = sum(1 for x in res if x[0]["interval"]["low"] <= d <= x[0]["interval"]["high"])
    ok = plo <= rej <= phi and clo <= cov <= chi
    lines.append("planted difference 0.5 sd of the differences, n=30: power (noncentral t, scipy) %.6f, rejected %d of %d (interval %d..%d); the 95%% interval covered the planted difference %d times (%d..%d)%s" % (power, rej, reps, plo, phi, cov, clo, chi, "" if ok else "  <-- FAILS"))
    if not ok:
        bad.append("power or coverage of the planted difference")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--kind", default="t2", choices=["t2", "paired"])
    ap.add_argument("--replicates", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    bad = []
    lines = []
    if a.kind == "t2":
        refused = {"pooled-for-welch": False, "df-plus-one": False, "one-sided": False, "normal": False}
    else:
        refused = {"independent-samples": False, "df-plus-one": False, "one-sided": False, "normal": False}
    try:
        (calibrate_t2 if a.kind == "t2" else calibrate_paired)(a, lines, bad, refused)
    except Exception as e:
        print("gate_test2_calibration: harness failure: %s" % e)
        return 2
    for name, was in refused.items():
        lines.append("control %s: %s" % (name, "refused by the calibration test" if was else "NOT refused"))
        if not was:
            bad.append("control %s was not refused: the gate cannot tell it from the program" % name)
    print("gate_test2_calibration %s (%d replicates per configuration, seed %d):" % (a.kind, a.replicates, a.seed))
    for ln in lines:
        print("  " + ln)
    print("gate_test2_calibration: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
