#!/usr/bin/env python3
"""S3.3 gates T-E (calibration) of docs/test.md sections 5 and 10 for `stats test chisq`:

    python3 scripts/gate_test3_calibration.py [--bin build/stats] [--replicates 4000] [--seed 1] [--workers 8]

Under independence (A and G drawn independently from fixed margins) and under a true fit (A drawn from the proportions), with every expected count at least 5, the rate at which `p < 0.05` lies inside the
two-sided 99.9% binomial interval for the number of replicates; and where the table is large enough for the p-values not to sit on a lattice (N of 2,000 or more), the p-values are uniform too (Kolmogorov-Smirnov,
p > 0.001). A table with an expected count under 5 is printed, not gated (the approximation is the thing known to fail there). CONTROLS read from the same replicates must be refused by the same test at some
configuration, so the gate is able to fail: the degrees of freedom r c - 1 (k for a fit) in place of (r - 1)(c - 1) (k - 1), and for tables of 2 x 2 Yates's correction.

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
from gate_test_calibration import ALPHA_TEST, interval  # noqa: E402

# (kind, N, margins of A, margins of G or the proportions text): independence tables have every expected count at least 5; the last of each kind is large enough for the KS test
CONFIGS = [
    ("ind", 40, [0.5, 0.5], [0.5, 0.5]),
    ("ind", 100, [0.5, 0.5], [0.5, 0.5]),
    ("ind", 300, [0.3, 0.3, 0.4], [0.2, 0.3, 0.5]),
    ("ind", 200, [0.5, 0.5], [0.2, 0.2, 0.2, 0.2, 0.2]),
    ("ind", 2000, [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]),
    ("ind", 3000, [0.5, 0.5], [0.4, 0.6]),
    ("fit", 100, [0.2, 0.3, 0.5], None),
    ("fit", 400, [0.1, 0.2, 0.3, 0.4], None),
    ("fit", 2000, [0.25, 0.25, 0.5], None),
]


def draw(r, margins):
    u = (r.next() % 1000000) / 1000000.0
    acc = 0.0
    for i, m in enumerate(margins):
        acc += m
        if u < acc:
            return i
    return len(margins) - 1


def one(binary, d, job):
    kind, a, g, props = job
    if kind == "ind":
        (d / "f.csv").write_text("A,G\n" + "\n".join("a%d,g%d" % (x, y) for x, y in zip(a, g)) + "\n")
        args = ["chisq", "--cols", "A", "--vs", "G", "f.csv"]
    else:
        (d / "f.csv").write_text("A\n" + "\n".join("a%d" % x for x in a) + "\n")
        args = ["chisq", "--cols", "A", "--expected", props, "f.csv"]
    p = subprocess.run([binary, "test", "--root", str(d)] + args, capture_output=True, text=True, timeout=120)
    if p.returncode != 0:
        raise RuntimeError("exit %d: %s" % (p.returncode, p.stdout[:600]))
    return json.loads(p.stdout)["data"]["tests"][0]


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


def judge(ps, lo, hi, with_ks):
    rej = sum(1 for p in ps if p < 0.05)
    ks = ss.kstest(ps, "uniform").pvalue
    return lo <= rej <= hi and (not with_ks or ks > ALPHA_TEST), rej, ks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--replicates", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    reps = a.replicates
    lo, hi = interval(reps, 0.05)
    bad, lines = [], []
    refused = {"df-rc-minus-one": False, "df-k": False, "yates": False}
    try:
        for ci, (kind, n, ma, mg) in enumerate(CONFIGS):
            r = SM64(a.seed * 1000 + ci)
            if kind == "ind":
                jobs = [(kind, [draw(r, ma) for _ in range(n)], [draw(r, mg) for _ in range(n)], None) for _ in range(reps)]
            else:
                props = ",".join("%g" % m for m in ma)
                jobs = [(kind, [draw(r, ma) for _ in range(n)], None, props) for _ in range(reps)]
            res = run_many(a.bin, jobs, a.workers)
            ps = [x["p_value"] for x in res if x["p_value"] is not None]
            below = [x for x in res if x["assumptions"].get("expected_below_5", 0) > 0]
            # gated when the DESIGN has every expected count at least 5 (a sample table may have a cell under 5 by chance: the rate of those is printed)
            design_min = n * min(ma) * (min(mg) if kind == "ind" else 1.0)
            gated = design_min >= 5
            with_ks = n >= 2000
            ok, rej, ks = judge(ps, lo, hi, with_ks)
            shape = "%d x %d" % (len(ma), len(mg)) if kind == "ind" else "k = %d" % len(ma)
            lines.append("%s %s N=%d: rejected %d of %d (interval %d..%d), KS p %.3f%s%s%s" % (kind, shape, n, rej, len(ps), lo, hi, ks, "" if with_ks else " (not gated: a lattice)", "" if ok else ("  <-- FAILS" if gated else "  (printed, not gated: the design has an expected count under 5)"), "; %d tables with a cell under 5" % len(below) if below else ""))
            if not ok and gated:
                bad.append("%s %s N=%d null calibration" % (kind, shape, n))
            xs = [x["statistic"] for x in res if x["statistic"] is not None]
            dfs = [x["df"] for x in res if x["statistic"] is not None]
            if kind == "ind":
                cp = [float(ss.chi2.sf(x, (len(ma)) * (len(mg)) - 1)) for x in xs]
            else:
                cp = [float(ss.chi2.sf(x, len(ma))) for x in xs]
            okc, _, _ = judge(cp, lo, hi, with_ks)
            if not okc:
                refused["df-rc-minus-one" if kind == "ind" else "df-k"] = True
            if kind == "ind" and len(ma) == 2 and len(mg) == 2:
                # Yates's correction on the same tables: the statistic is not in the answer, so recompute it from the counts of the replicate
                yp = []
                for job in jobs:
                    _, av, gv, _ = job
                    t = [[0, 0], [0, 0]]
                    for x, y in zip(av, gv):
                        t[x][y] += 1
                    yp.append(float(ss.chi2_contingency(t, correction=True)[1]))
                oky, _, _ = judge(yp, lo, hi, with_ks)
                if not oky:
                    refused["yates"] = True
    except Exception as e:
        print("gate_test3_calibration: harness failure: %s" % e)
        return 2
    for name, was in refused.items():
        lines.append("control %s: %s" % (name, "refused by the calibration test" if was else "NOT refused"))
        if not was:
            bad.append("control %s was not refused: the gate cannot tell it from the program" % name)
    print("gate_test3_calibration chisq (%d replicates per configuration, seed %d):" % (reps, a.seed))
    for ln in lines:
        print("  " + ln)
    print("gate_test3_calibration: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
