#!/usr/bin/env python3
"""S3 gates T-E (calibration), T-F (power and coverage) of docs/test.md section 5, kind by kind (`t1` so far).

    python3 scripts/gate_test_calibration.py [--bin build/stats] [--kind t1] [--replicates 4000] [--seed 1] [--workers 8]

t1. Under a true null (n normal values of mean 0, seeded) and for n = 3, 5, 10, 30, 100, 1,000 the rate at which `p < 0.05` lies inside the two-sided 99.9% binomial interval for the number of replicates, and
the p-values are uniform (Kolmogorov-Smirnov, p > 0.001). The 95% interval (`--level 0.95`) covers the true mean at the same rate, tested the same way, and the 90% interval likewise: that exercises the
critical value and the standard error, which the p-value alone does not. CONTROLS from the same replicates (a one-sided p-value, the normal in place of Student's t, `n` degrees of freedom) must be refused
by the same test at some n, so the gate is able to fail. Power: with a planted mean of 0.5 standard deviations and n = 30 the rate at which `p < 0.05` equals the power of the noncentral t (scipy, an
independent reference), inside the binomial interval; the coverage of the planted mean is again nominal.

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

NS = (3, 5, 10, 30, 100, 1000)
ALPHA_TEST = 0.001


def normal(r):
    u1 = max(r.u(), 1e-300)
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * r.u())


def one(binary, d, xs):
    f = d / "f.csv"
    f.write_text("x\n" + "\n".join(repr(v) for v in xs) + "\n")
    out = []
    for level in ("0.95", "0.9"):
        p = subprocess.run([binary, "test", "t1", "--root", str(d), "--cols", "x", "--mu0", "0", "--level", level, "f.csv"], capture_output=True, text=True, timeout=120)
        if p.returncode != 0:
            raise RuntimeError("exit %d: %s" % (p.returncode, p.stdout[:200]))
        t = json.loads(p.stdout)["data"]["tests"][0]
        out.append(t)
    return out


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


def interval(reps, p):
    return int(ss.binom.ppf(ALPHA_TEST / 2, reps, p)), int(ss.binom.ppf(1 - ALPHA_TEST / 2, reps, p))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--kind", default="t1", choices=["t1"])
    ap.add_argument("--replicates", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    bad = []
    lines = []
    refused = {"one-sided": False, "normal": False, "df-n": False}
    reps = a.replicates
    lo, hi = interval(reps, 0.05)
    clo, chi = interval(reps, 0.95)
    c9lo, c9hi = interval(reps, 0.9)
    try:
        for n in NS:
            r = SM64(a.seed * 1000 + n)
            jobs = [[normal(r) for _ in range(n)] for _ in range(reps)]
            res = run_many(a.bin, jobs, a.workers)
            ps = [t95["p_value"] for t95, _ in res]
            ts = [t95["statistic"] for t95, _ in res]
            rej = sum(1 for p in ps if p < 0.05)
            ks = ss.kstest(ps, "uniform").pvalue
            cov95 = sum(1 for t95, _ in res if t95["interval"]["low"] <= 0.0 <= t95["interval"]["high"])
            cov90 = sum(1 for _, t90 in res if t90["interval"]["low"] <= 0.0 <= t90["interval"]["high"])
            ok = lo <= rej <= hi and ks > ALPHA_TEST and clo <= cov95 <= chi and c9lo <= cov90 <= c9hi
            lines.append("n=%d: rejected %d of %d (interval %d..%d), KS p %.3f; the 95%% interval covered %d (%d..%d), the 90%% interval %d (%d..%d)%s" % (n, rej, reps, lo, hi, ks, cov95, clo, chi, cov90, c9lo, c9hi, "" if ok else "  <-- FAILS"))
            if not ok:
                bad.append("n=%d null calibration or coverage" % n)
            for name in refused:
                cp = []
                for p, t in zip(ps, ts):
                    if name == "one-sided":
                        cp.append(p / 2.0)
                    elif name == "normal":
                        cp.append(2.0 * ss.norm.sf(abs(t)))
                    else:
                        cp.append(2.0 * ss.t.sf(abs(t), n))
                crej = sum(1 for p in cp if p < 0.05)
                cks = ss.kstest(cp, "uniform").pvalue
                if not (lo <= crej <= hi and cks > ALPHA_TEST):
                    refused[name] = True
        # power and the coverage of a planted mean
        n, delta = 30, 0.5
        pr = SM64(a.seed * 7919)
        jobs = [[delta + normal(pr) for _ in range(n)] for _ in range(reps)]
        res = run_many(a.bin, jobs, a.workers)
    except Exception as e:
        print("gate_test_calibration: harness failure: %s" % e)
        return 2
    tc = ss.t.ppf(0.975, n - 1)
    nc = delta * math.sqrt(n)
    power = float(ss.nct.sf(tc, n - 1, nc) + ss.nct.cdf(-tc, n - 1, nc))
    rej = sum(1 for t95, _ in res if t95["p_value"] < 0.05)
    plo, phi = interval(reps, power)
    cov = sum(1 for t95, _ in res if t95["interval"]["low"] <= delta <= t95["interval"]["high"])
    ok = plo <= rej <= phi and clo <= cov <= chi
    lines.append("planted mean 0.5 sd, n=30: power (noncentral t, scipy) %.6f, rejected %d of %d (interval %d..%d); the 95%% interval covered the planted mean %d times (%d..%d)%s" % (power, rej, reps, plo, phi, cov, clo, chi, "" if ok else "  <-- FAILS"))
    if not ok:
        bad.append("power or coverage of the planted mean")
    for name, was in refused.items():
        lines.append("control %s: %s" % (name, "refused by the calibration test" if was else "NOT refused"))
        if not was:
            bad.append("control %s was not refused: the gate cannot tell it from the program" % name)
    print("gate_test_calibration t1 (%d null replicates per n, seed %d):" % (reps, a.seed))
    for ln in lines:
        print("  " + ln)
    print("gate_test_calibration: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
