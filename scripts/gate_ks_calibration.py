#!/usr/bin/env python3
"""S4.5 gate K-D (calibration) and K-E (power, printed) of docs/rank.md section 7 for `stats test ks`:

    python3 scripts/gate_ks_calibration.py [--bin build/stats] [--replicates 4000] [--seed 1] [--workers 8]

Under a true null (two groups of one distribution, seeded: normal values, so no two values equal) for sizes from 5 to 300, the rate at which `p < 0.05` lies inside the two-sided 99.9% binomial interval for the number of replicates
around the EXACT SIZE OF THE RULE THAT WAS USED: the proportion of the arrangements whose distance d is at least the smallest d the rule rejects (the exact p-value, or the asymptotic Q_KS(lambda) at the exact lambda, both
monotone in d), counted by the oracle's recurrence, which is exact at every size. Neither rule has size 0.05: the exact test is conservative on a lattice of attainable p-values and the asymptotic one is conservative at small
sizes (no Stephens correction); the sizes are printed beside 0.05. With ties (a lattice of 5 values, a binary one) the conditional distribution is not the one counted, and the rates are printed, not gated. CONTROLS from the same
replicates must be refused by the same test at some configuration, so the gate is able to fail: 0.05 as the size of either rule, and the one-sided p-value. Power against a planted shift: printed.

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

import mpmath

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
from gate_test_calibration import interval, normal  # noqa: E402
import oracle_test as ot  # noqa: E402

CONFIGS = [
    ((5, 5), "normal"),
    ((8, 12), "normal"),
    ((10, 10), "normal"),
    ((20, 20), "normal"),
    ((30, 30), "normal"),
    ((31, 30), "normal"),
    ((100, 100), "normal"),
    ((100, 200), "normal"),
    ((300, 300), "normal"),
    ((20, 20), "lattice"),
    ((30, 30), "binary"),
]


def write(d, groups):
    rows = ["x,g"]
    for gi, vals in enumerate(groups):
        rows += ["%r,g%d" % (v, gi) for v in vals]
    (d / "f.csv").write_text("\n".join(rows) + "\n")


def one(binary, d, groups):
    write(d, groups)
    out = []
    for extra in ([], ["--exact", "no"]):
        p = subprocess.run([binary, "test", "ks", "--root", str(d), "--cols", "x", "--group", "g"] + extra + ["f.csv"], capture_output=True, text=True, timeout=120)
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


def q_ks(lam2):
    x = mpmath.sqrt(mpmath.mpf(lam2.numerator) / lam2.denominator)
    with mpmath.workdps(120):
        tot = mpmath.mpf(0)
        j = 1
        while True:
            term = mpmath.exp(-2 * j * j * x * x)
            tot += term if j % 2 else -term
            if term < mpmath.mpf(10) ** -130:
                break
            j += 1
        return min(mpmath.mpf(1), 2 * tot)


def size_of_rule(na, nb, exact):
    """the probability, under the null with no two values equal, that the rule's p-value is below 0.05: the proportion of arrangements with d at least the smallest rejected d"""
    nn = na + nb
    total = math.comb(nn, na)
    for d in range(1, na * nb + 1):
        if exact:
            p = Fraction(ot.ks_exceed(na, nb, d), total)
            hit = p < Fraction(5, 100)
        else:
            hit = q_ks(Fraction(d * d, nn * na * nb)) < mpmath.mpf("0.05")
        if hit:
            return float(Fraction(ot.ks_exceed(na, nb, d), total)), d
    return 0.0, na * nb + 1


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
    refused = {"nominal-size": False, "one-sided": False}
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
            for pos, label in ((0, "default"), (1, "asymptotic")):
                methods = set(x[pos]["p_method"] for x in res)
                ps = [x[pos]["p_value"] for x in res if x[pos]["p_value"] is not None]
                rej = sum(1 for p in ps if p < 0.05)
                if kind == "normal":
                    exact = methods == {"exact"}
                    size, dstar = size_of_rule(na, nb, exact)
                    ok, lo, hi = inside(rej, len(ps), size)
                    lines.append("n=%s %s %s (%s): exact size of the rule %.5f (d >= %d), rejected %d of %d (interval %d..%d)%s" % (sizes, kind, label, "exact" if exact else "asymptotic", size, dstar, rej, len(ps), lo, hi, "" if ok else "  <-- FAILS"))
                    if not ok:
                        bad.append("%s n=%s calibration" % (label, sizes))
                    okf, _, _ = inside(rej, len(ps), 0.05)
                    if not okf:
                        refused["nominal-size"] = True
                    oks, _, _ = inside(sum(1 for p in ps if p / 2 < 0.05), len(ps), size)
                    if not oks:
                        refused["one-sided"] = True
                else:
                    lines.append("n=%s %s %s (%s): rejected %d of %d (printed, not gated: ties change the conditional distribution)" % (sizes, kind, label, "/".join(sorted(methods)), rej, len(ps)))
        pr = SM64(a.seed * 7919)
        jobs = [[[0.5 + normal(pr) for _ in range(20)], [normal(pr) for _ in range(30)]] for _ in range(reps)]
        res = run_many(a.bin, jobs, a.workers)
        lines.append("planted shift 0.5 sd, n=(20,30): rejected %d of %d (printed, not gated)" % (sum(1 for x in res if x[0]["p_value"] < 0.05), reps))
    except Exception as e:
        print("gate_ks_calibration: harness failure: %s" % e)
        return 2
    for name, was in refused.items():
        lines.append("control %s: %s" % (name, "refused by the calibration test" if was else "NOT refused"))
        if not was:
            bad.append("control %s was not refused: the gate cannot tell it from the program" % name)
    print("gate_ks_calibration (%d replicates per configuration, seed %d):" % (reps, a.seed))
    for ln in lines:
        print("  " + ln)
    print("gate_ks_calibration: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
