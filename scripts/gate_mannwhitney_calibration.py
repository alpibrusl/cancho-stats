#!/usr/bin/env python3
"""S4.3 gate R-D (calibration) and R-E (power, printed) of docs/rank.md section 5 for `stats test mannwhitney`:

    python3 scripts/gate_mannwhitney_calibration.py [--bin build/stats] [--replicates 4000] [--seed 1] [--workers 8]

Under a true null (two groups of one distribution, seeded: normal values, a lattice of 5 values and a binary one for heavy ties), for sizes from 5 to 1,000, the rate at which `p < 0.05` lies inside the two-sided 99.9% binomial
interval for the number of replicates. WHERE THE P-VALUE IS EXACT (no ties, at most 60 values) the interval is around the EXACT SIZE OF THE TEST AS DEFINED, the sum of the probabilities of the values of U whose
two-sided p is below 0.05, computed by the oracle's counts (the exact test is conservative: its size is below 0.05, said and printed); where it is the normal approximation the interval is around 0.05 and the configurations
are those with 10 or more values in both groups (the others are printed). The p-values are uniform (Kolmogorov-Smirnov) only for the normal approximation without ties. CONTROLS from the same replicates must be refused by
the same test at some configuration, so the gate is able to fail: the normal p-value where the exact one is the answer (compared with the exact size), the normal p-value without the continuity correction, the one-sided
p-value, and 0.05 as the size of an exact test. Power against a planted shift: printed.

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
from fractions import Fraction

import scipy.stats as ss

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
from gate_test_calibration import ALPHA_TEST, interval, normal  # noqa: E402
import oracle_test as ot  # noqa: E402

CONFIGS = [
    ((5, 5), "normal"),
    ((8, 12), "normal"),
    ((20, 20), "normal"),
    ((30, 30), "normal"),
    ((31, 30), "normal"),
    ((100, 100), "normal"),
    ((300, 700), "normal"),
    ((10, 50), "normal"),
    ((20, 20), "lattice"),
    ((30, 30), "binary"),
    ((100, 200), "binary"),
]


def write(d, groups):
    rows = ["x,g"]
    for gi, vals in enumerate(groups):
        rows += ["%r,g%d" % (v, gi) for v in vals]
    (d / "f.csv").write_text("\n".join(rows) + "\n")


def one(binary, d, groups):
    write(d, groups)
    out = []
    for extra in ([], ["--exact", "no"], ["--exact", "no", "--continuity", "none"]):
        p = subprocess.run([binary, "test", "mannwhitney", "--root", str(d), "--cols", "x", "--group", "g"] + extra + ["f.csv"], capture_output=True, text=True, timeout=120)
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


def exact_size(na, nb):
    """the probability, under the null, that the exact two-sided p is below 0.05: the sum of the probabilities of the values of U whose p is"""
    counts = ot.mw_counts(na, nb)
    total = sum(counts)
    size = Fraction(0)
    for u in range(len(counts)):
        le = sum(counts[:u + 1])
        ge = sum(counts[u:])
        p = min(Fraction(1), Fraction(2 * min(le, ge), total))
        if float(p) < 0.05:
            size += Fraction(counts[u], total)
    return float(size)


def inside(rej, reps, size):
    lo, hi = interval(reps, size)
    return lo <= rej <= hi, lo, hi


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
    refused = {"normal-for-exact": False, "no-continuity-correction": False, "one-sided": False, "nominal-size-for-exact": False}
    try:
        for ci, (sizes, kind) in enumerate(CONFIGS):
            r = SM64(a.seed * 1000 + ci)
            na, nb = sizes
            if kind == "normal":
                jobs = [[[normal(r) for _ in range(n)] for n in sizes] for _ in range(reps)]
            else:
                top = 5 if kind == "lattice" else 2
                jobs = [[[float(int(r.next() % top)) for _ in range(n)] for n in sizes] for _ in range(reps)]
            res = run_many(a.bin, jobs, a.workers)
            methods = set(x[0]["p_method"] for x in res)
            ps = [x[0]["p_value"] for x in res if x[0]["p_value"] is not None]
            rej = sum(1 for p in ps if p < 0.05)
            if methods == {"exact"}:
                size = exact_size(na, nb)
                ok, lo, hi = inside(rej, len(ps), size)
                lines.append("n=%s %s exact: exact size of the test %.5f, rejected %d of %d (interval %d..%d)%s" % (sizes, kind, size, rej, len(ps), lo, hi, "" if ok else "  <-- FAILS"))
                if not ok:
                    bad.append("exact n=%s calibration" % (sizes,))
                # controls where they differ from the exact size
                nrm = [x[1]["p_value"] for x in res]
                okn, _, _ = inside(sum(1 for p in nrm if p < 0.05), reps, size)
                if not okn:
                    refused["normal-for-exact"] = True
                okf, _, _ = inside(rej, reps, 0.05)
                if not okf:
                    refused["nominal-size-for-exact"] = True
                oks, _, _ = inside(sum(1 for p in ps if p / 2 < 0.05), len(ps), size)
                if not oks:
                    refused["one-sided"] = True
            else:
                lo, hi = interval(reps, 0.05)
                gated = min(sizes) >= 10
                ks = ss.kstest(ps, "uniform").pvalue if kind == "normal" else 1.0
                ok = lo <= rej <= hi and ks > ALPHA_TEST
                lines.append("n=%s %s normal (tie-corrected, continuity-corrected): rejected %d of %d (interval %d..%d)%s%s" % (sizes, kind, rej, len(ps), lo, hi, (", KS p %.3f" % ks) if kind == "normal" else "", ("" if ok else "  <-- FAILS") if gated else ("" if ok else "  (printed, not gated: a group below 10)")))
                if not ok and gated:
                    bad.append("normal n=%s %s calibration" % (sizes, kind))
                if gated:
                    nc = [x[2]["p_value"] for x in res]
                    okc, _, _ = inside(sum(1 for p in nc if p < 0.05), reps, 0.05)
                    if not okc or ss.kstest(nc, "uniform").pvalue <= ALPHA_TEST and kind == "normal":
                        refused["no-continuity-correction"] = True
                    oko, _, _ = inside(sum(1 for p in ps if p / 2 < 0.05), len(ps), 0.05)
                    if not oko:
                        refused["one-sided"] = True
        # power: printed
        pr = SM64(a.seed * 7919)
        jobs = [[[0.5 + normal(pr) for _ in range(20)], [normal(pr) for _ in range(30)]] for _ in range(reps)]
        res = run_many(a.bin, jobs, a.workers)
        lines.append("planted shift 0.5 sd, n=(20,30): rejected %d of %d (printed, not gated)" % (sum(1 for x in res if x[0]["p_value"] < 0.05), reps))
    except Exception as e:
        print("gate_mannwhitney_calibration: harness failure: %s" % e)
        return 2
    for name, was in refused.items():
        lines.append("control %s: %s" % (name, "refused by the calibration test" if was else "NOT refused"))
        if not was:
            bad.append("control %s was not refused: the gate cannot tell it from the program" % name)
    print("gate_mannwhitney_calibration (%d replicates per configuration, seed %d):" % (reps, a.seed))
    for ln in lines:
        print("  " + ln)
    print("gate_mannwhitney_calibration: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
