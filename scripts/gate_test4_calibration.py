#!/usr/bin/env python3
"""S3.4 gates T-E (calibration) and T-F (coverage) of docs/test.md sections 5 and 11 for `stats test prop`:

    python3 scripts/gate_test4_calibration.py [--bin build/stats] [--replicates 4000] [--seed 1] [--workers 8]

The counts of a proportion are discrete, so the size of the z test is not the nominal 5% (at n = 50 and p0 = 1/2 it is 6.6%) and the coverage of a score interval oscillates around its level. The gate therefore compares the
rates the program shows with the EXACT probabilities of its own definitions, computed by summing the binomial probabilities over every outcome (numpy, scipy: independent of the program): for the one-sample z at p0 = the true
proportion and for the pooled two-sample z at p1 = p2, the rate at which `p < 0.05` lies inside the two-sided 99.9% binomial interval around the exact size; and the coverage of the Wilson interval (of the true proportion) and of
Newcombe's (of 0 for equal proportions) lies inside the interval around the exact coverage. The nominal values are printed beside the exact ones. CONTROLS must be refused by the same test at some configuration, so the gate is
able to fail: the one-sided p-value and Yates's continuity correction (both against the exact size of the uncorrected test), and Wald's interval (its coverage against the exact coverage of Wilson's).

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

import numpy as np
import scipy.stats as ss

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
from gate_test_calibration import interval  # noqa: E402

ONE = [(50, 0.5), (100, 0.3), (400, 0.1), (1000, 0.02), (2000, 0.5)]
TWO = [(100, 100, 0.3), (50, 200, 0.5), (400, 100, 0.1), (1000, 1000, 0.05)]
Q = float(ss.norm.isf(0.025))


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


def call(binary, d, args):
    p = subprocess.run([binary, "test", "prop", "--root", str(d)] + args, capture_output=True, text=True, timeout=120)
    if p.returncode != 0:
        raise RuntimeError("exit %d: %s" % (p.returncode, p.stdout[:400]))
    return json.loads(p.stdout)["data"]["tests"][0]


def one_job(binary, d, job):
    kind, xs, p0 = job
    if kind == "one":
        (d / "f.csv").write_text("y\n" + "\n".join("s" if v else "f" for v in xs) + "\n")
        return call(binary, d, ["--cols", "y", "--success", "s", "--p0", "%g" % p0, "f.csv"])
    a, b = xs
    rows = ["s,a" if v else "f,a" for v in a] + ["s,b" if v else "f,b" for v in b]
    (d / "f.csv").write_text("y,g\n" + "\n".join(rows) + "\n")
    return call(binary, d, ["--cols", "y", "--success", "s", "--group", "g", "f.csv"])


def wilson(x, n, q):
    x = np.asarray(x, dtype=float)
    c = (x + q * q / 2) / (n + q * q)
    h = q * np.sqrt(x * (n - x) / n + q * q / 4) / (n + q * q)
    return np.maximum(c - h, 0.0), np.minimum(c + h, 1.0)


def exact_one(n, p0, ctl=False):
    """(size of the z test at p0, coverage of Wilson's interval of p0, coverage of Wald's), by summing the binomial probabilities"""
    xs = np.arange(n + 1)
    pm = ss.binom.pmf(xs, n, p0)
    z = (xs - n * p0) / math.sqrt(n * p0 * (1 - p0))
    size = float(pm[2 * ss.norm.sf(np.abs(z)) < 0.05].sum())
    lo, hi = wilson(xs, n, Q)
    cov = float(pm[(lo <= p0) & (p0 <= hi)].sum())
    ph = xs / n
    half = Q * np.sqrt(ph * (1 - ph) / n)
    wald = float(pm[(ph - half <= p0) & (p0 <= ph + half)].sum())
    yates = float(pm[2 * ss.norm.sf(np.maximum(np.abs(xs - n * p0) - 0.5, 0) / math.sqrt(n * p0 * (1 - p0))) < 0.05].sum())
    return size, cov, wald, yates


def exact_two(n1, n2, p):
    x1 = np.arange(n1 + 1)[:, None]
    x2 = np.arange(n2 + 1)[None, :]
    pm = ss.binom.pmf(x1, n1, p) * ss.binom.pmf(x2, n2, p)
    nn, xx = n1 + n2, x1 + x2
    ok = (xx > 0) & (xx < nn)
    pp = np.where(ok, xx / nn, 0.5)
    d = x1 / n1 - x2 / n2
    z = d / np.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
    rej = ok & (2 * ss.norm.sf(np.abs(z)) < 0.05)
    l1, u1 = wilson(x1, n1, Q)
    l2, u2 = wilson(x2, n2, Q)
    p1, p2 = x1 / n1, x2 / n2
    lo = d - np.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)
    hi = d + np.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
    cov = (lo <= 0) & (0 <= hi)
    corr = 0.5 * (1 / n1 + 1 / n2)
    zc = np.maximum(np.abs(d) - corr, 0) / np.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
    yates = ok & (2 * ss.norm.sf(zc) < 0.05)
    return float(pm[rej].sum()), float(pm[cov].sum()), float(pm[yates].sum())


def inside(count, reps, prob):
    lo, hi = interval(reps, prob)
    return lo <= count <= hi, "%d (%d..%d)" % (count, lo, hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--replicates", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    reps = a.replicates
    bad, lines = [], []
    refused = {"one-sided": False, "yates": False, "wald": False}
    try:
        for ci, (n, p0) in enumerate(ONE):
            r = SM64(a.seed * 1000 + ci)
            jobs = [("one", [(r.next() % 1000000) / 1000000.0 < p0 for _ in range(n)], p0) for _ in range(reps)]
            res = run_many(one_job, a.bin, jobs, a.workers)
            size, cov, wald, yates = exact_one(n, p0)
            rej = sum(1 for x in res if x["p_value"] is not None and x["p_value"] < 0.05)
            covered = sum(1 for x, j in zip(res, jobs) if x["interval"]["low"] <= p0 <= x["interval"]["high"])
            ok1, s1 = inside(rej, reps, size)
            ok2, s2 = inside(covered, reps, cov)
            lines.append("one sample n=%d p0=%g: rejected %s, exact size %.4f (nominal 0.05); Wilson covered %s, exact coverage %.4f (nominal 0.95)%s" % (n, p0, s1, size, s2, cov, "" if ok1 and ok2 else "  <-- FAILS"))
            if not (ok1 and ok2):
                bad.append("one sample n=%d p0=%g" % (n, p0))
            # controls: the same replicates read wrongly
            if not inside(sum(1 for x in res if x["p_value"] is not None and x["p_value"] / 2 < 0.05), reps, size)[0]:
                refused["one-sided"] = True
            yrej = sum(1 for x, j in zip(res, jobs) if abs(sum(j[1]) - n * p0) - 0.5 > 0 and 2 * ss.norm.sf((abs(sum(j[1]) - n * p0) - 0.5) / math.sqrt(n * p0 * (1 - p0))) < 0.05)
            if not inside(yrej, reps, size)[0]:
                refused["yates"] = True
            wcov = 0
            for j in jobs:
                x = sum(j[1])
                ph = x / n
                half = Q * math.sqrt(ph * (1 - ph) / n)
                wcov += 1 if ph - half <= p0 <= ph + half else 0
            if not inside(wcov, reps, cov)[0]:
                refused["wald"] = True
        for ci, (n1, n2, p) in enumerate(TWO):
            r = SM64(a.seed * 2000 + ci)
            jobs = [("two", ([(r.next() % 1000000) / 1000000.0 < p for _ in range(n1)], [(r.next() % 1000000) / 1000000.0 < p for _ in range(n2)]), None) for _ in range(reps)]
            res = run_many(one_job, a.bin, jobs, a.workers)
            size, cov, yates = exact_two(n1, n2, p)
            rej = sum(1 for x in res if x["p_value"] is not None and x["p_value"] < 0.05)
            covered = sum(1 for x in res if x["interval"] is not None and x["interval"]["low"] <= 0 <= x["interval"]["high"])
            ok1, s1 = inside(rej, reps, size)
            ok2, s2 = inside(covered, reps, cov)
            lines.append("two samples n=(%d,%d) p=%g: rejected %s, exact size %.4f; Newcombe covered 0 %s, exact coverage %.4f%s" % (n1, n2, p, s1, size, s2, cov, "" if ok1 and ok2 else "  <-- FAILS"))
            if not (ok1 and ok2):
                bad.append("two samples n=(%d,%d) p=%g" % (n1, n2, p))
            if not inside(sum(1 for x in res if x["p_value"] is not None and x["p_value"] / 2 < 0.05), reps, size)[0]:
                refused["one-sided"] = True
            if not inside(rej, reps, yates)[0]:
                refused["yates"] = True
    except Exception as e:
        print("gate_test4_calibration: harness failure: %s" % e)
        return 2
    for name, was in refused.items():
        lines.append("control %s: %s" % (name, "refused by the calibration test" if was else "NOT refused"))
        if not was:
            bad.append("control %s was not refused: the gate cannot tell it from the program" % name)
    print("gate_test4_calibration prop (%d replicates per configuration, seed %d):" % (reps, a.seed))
    for ln in lines:
        print("  " + ln)
    print("gate_test4_calibration: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
