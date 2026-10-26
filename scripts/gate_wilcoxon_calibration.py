#!/usr/bin/env python3
"""S4.4 gate R-D (calibration) and R-E (power, printed) of docs/rank.md section 5 for `stats test wilcoxon`:

    python3 scripts/gate_wilcoxon_calibration.py [--bin build/stats] [--replicates 4000] [--seed 1] [--workers 8]

Under a true null (n values from a distribution symmetric about 0, seeded: normal, and a lattice of the integers -2 .. 2 for zeros and heavy ties), for n from 5 to 300, the rate at which `p < 0.05` of the one-sample test
(`--cols x --mu0 0`) lies inside the two-sided 99.9% binomial interval for the number of replicates. WHERE THE P-VALUE IS EXACT (no ties, no zeros, at most 50 values) the interval is around the EXACT SIZE OF THE TEST AS
DEFINED, the sum of the probabilities of the values of W+ whose two-sided p is below 0.05, computed from the oracle's counts (conservative; printed); where it is the normal approximation the interval is around 0.05 and
the configurations are those with 20 or more nonzero values (the others are printed). On the lattice the default method for zeros (`wilcox`, conditional on the nonzero differences) and the other two are all judged against
0.05. The p-values are uniform (Kolmogorov-Smirnov) only for the normal approximation without ties. CONTROLS from the same replicates must be refused at some configuration, so the gate is able to fail: the normal p-value
WITHOUT the continuity correction where the exact one is the answer (compared with the exact size), the one-sided p-value, 0.05 as the size of an exact test. Power against a planted shift: printed.

Exit 0: every test passes and every control is refused. Exit 1: a test fails. Exit 2: the harness failed.
"""
import argparse
import concurrent.futures
import json
import pathlib
import subprocess
import sys
import tempfile
from fractions import Fraction

import scipy.stats as ss

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
from gate_test_calibration import ALPHA_TEST, interval, normal  # noqa: E402
import oracle_test as ot  # noqa: E402

CONFIGS = [(5, "normal"), (8, "normal"), (12, "normal"), (20, "normal"), (35, "normal"), (50, "normal"), (51, "normal"), (100, "normal"), (300, "normal"), (40, "lattice"), (100, "lattice")]


def write(d, xs):
    (d / "f.csv").write_text("x\n" + "\n".join("%r" % v for v in xs) + "\n")


def one(binary, d, xs):
    write(d, xs)
    out = []
    for extra in ([], ["--exact", "no", "--continuity", "none"], ["--zeros", "pratt"], ["--zeros", "zsplit"]):
        p = subprocess.run([binary, "test", "wilcoxon", "--root", str(d), "--cols", "x", "--mu0", "0"] + extra + ["f.csv"], capture_output=True, text=True, timeout=120)
        if p.returncode != 0:
            raise RuntimeError("exit %d: %s" % (p.returncode, p.stdout[:200]))
        out.append(json.loads(p.stdout)["data"]["tests"][0])
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


def exact_size(n):
    counts = ot.wx_counts(n)
    total = 2 ** n
    size = Fraction(0)
    for w in range(len(counts)):
        le = sum(counts[:w + 1])
        ge = sum(counts[w:])
        p = min(Fraction(1), Fraction(2 * min(le, ge), total))
        if float(p) < 0.05:
            size += Fraction(counts[w], total)
    return float(size)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--replicates", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    reps = a.replicates
    bad = []
    lines = []
    refused = {"no-continuity-correction": False, "one-sided": False, "nominal-size-for-exact": False}
    try:
        for ci, (n, kind) in enumerate(CONFIGS):
            r = SM64(a.seed * 1000 + ci)
            if kind == "normal":
                jobs = [[normal(r) for _ in range(n)] for _ in range(reps)]
            else:
                jobs = [[float(int(r.next() % 5) - 2) for _ in range(n)] for _ in range(reps)]
            res = run_many(a.bin, jobs, a.workers)
            ps = [x[0]["p_value"] for x in res if x[0]["p_value"] is not None]
            rej = sum(1 for p in ps if p < 0.05)
            methods = set(x[0]["p_method"] for x in res)
            if methods == {"exact"}:
                size = exact_size(n)
                lo, hi = interval(len(ps), size)
                ok = lo <= rej <= hi
                lines.append("n=%d %s exact: exact size of the test %.5f, rejected %d of %d (interval %d..%d)%s" % (n, kind, size, rej, len(ps), lo, hi, "" if ok else "  <-- FAILS"))
                if not ok:
                    bad.append("exact n=%d calibration" % n)
                nc = [x[1]["p_value"] for x in res]
                lo2, hi2 = interval(reps, size)
                if not lo2 <= sum(1 for p in nc if p < 0.05) <= hi2:
                    refused["no-continuity-correction"] = True
                if not interval(reps, 0.05)[0] <= rej <= interval(reps, 0.05)[1]:
                    refused["nominal-size-for-exact"] = True
                if not lo <= sum(1 for p in ps if p / 2 < 0.05) <= hi:
                    refused["one-sided"] = True
            else:
                lo, hi = interval(len(ps), 0.05)
                gated = n >= 20
                ks = ss.kstest(ps, "uniform").pvalue if kind == "normal" else 1.0
                ok = lo <= rej <= hi and ks > ALPHA_TEST
                lines.append("n=%d %s normal (continuity-corrected%s): rejected %d of %d (interval %d..%d)%s%s" % (n, kind, ", tie-corrected" if kind == "lattice" else "", rej, len(ps), lo, hi, (", KS p %.3f" % ks) if kind == "normal" else "", ("" if ok else "  <-- FAILS") if gated else ("" if ok else "  (printed, not gated: fewer than 20)")))
                if not ok and gated:
                    bad.append("normal n=%d %s calibration" % (n, kind))
                if kind == "lattice":
                    for idx, nm in ((2, "pratt"), (3, "zsplit")):
                        pz = [x[idx]["p_value"] for x in res if x[idx]["p_value"] is not None]
                        rz = sum(1 for p in pz if p < 0.05)
                        loz, hiz = interval(len(pz), 0.05)
                        lines.append("   zeros %s: rejected %d of %d (interval %d..%d)" % (nm, rz, len(pz), loz, hiz))
                elif gated:
                    if not lo <= sum(1 for p in ps if p / 2 < 0.05) <= hi:
                        refused["one-sided"] = True
        pr = SM64(a.seed * 7919)
        jobs = [[0.5 + normal(pr) for _ in range(30)] for _ in range(reps)]
        res = run_many(a.bin, jobs, a.workers)
        lines.append("planted shift 0.5 sd, n=30: rejected %d of %d (printed, not gated)" % (sum(1 for x in res if x[0]["p_value"] < 0.05), reps))
    except Exception as e:
        print("gate_wilcoxon_calibration: harness failure: %s" % e)
        return 2
    for name, was in refused.items():
        lines.append("control %s: %s" % (name, "refused by the calibration test" if was else "NOT refused"))
        if not was:
            bad.append("control %s was not refused: the gate cannot tell it from the program" % name)
    print("gate_wilcoxon_calibration (%d replicates per configuration, seed %d):" % (reps, a.seed))
    for ln in lines:
        print("  " + ln)
    print("gate_wilcoxon_calibration: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
