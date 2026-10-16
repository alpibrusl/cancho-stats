#!/usr/bin/env python3
"""S1-A part 2 (docs/describe.md section 5): the exact power sums sum(x^k), k = 1 to 4, of lib/acc/psum.cho against `fractions.Fraction`, rounded once, on vectors of the families of measure_acc.py
(uniform, every finite exponent, huge, tiny, subnormal, cancellation, constant, edges of the range) and the extremes (the largest double, the smallest subnormal, alternating signs). Also the row
order: every vector is run again in a shuffled order and must print the same bytes.

    python3 scripts/gate_psum.py [--driver build/psum-check] [--vectors 300] [--seed 1]

Exit 0: every answer is the exactly rounded value (or code 1 exactly when the exact sum is beyond the largest double) and every shuffle agrees. Exit 1: a difference. Exit 2: the harness failed.
"""
import argparse
import math
import pathlib
import struct
import subprocess
import sys
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64, hx  # noqa: E402
import measure_acc as ma  # noqa: E402

MAXF = float.fromhex("0x1.fffffffffffffp+1023")
TINY = 5e-324


def expect(xs, k):
    s = sum((Fraction(x) ** k for x in xs), Fraction(0))
    try:
        return "%016x 0" % struct.unpack("<Q", struct.pack("<d", float(s)))[0]
    except OverflowError:
        return "%016x 1" % 0


def families(r, count):
    gens = [ma.fam_uniform, ma.mk(-1074, 1023), ma.mk(400, 1023), ma.mk(-1074, -480), ma.mk(-300, 300), ma.mk(-60, 60), ma.fam_cancel, ma.fam_const, ma.mk(-1074, -1023), ma.fam_edge]
    out = []
    for i in range(count):
        g = gens[i % len(gens)]
        n = [1, 2, 3, 7, 30, 200][i % 6]
        xs = g(r, n)[0]
        out.append(xs)
    out.append([MAXF, -MAXF, MAXF])
    out.append([TINY, -TINY, TINY, TINY])
    out.append([TINY] * 50)
    out.append([MAXF, MAXF])
    out.append([2.0 ** 1000, 2.0 ** 1000])
    out.append([float((-1) ** i) * 3.0 for i in range(100)])
    out.append([1e150, 1e-150, -1e150, 1e-150])
    out.append([MAXF] * 257)  # 257 rows put the sum just above 2^2106 units: a limb count without headroom wraps it to a finite double
    out.append([MAXF] * 4096)  # the headroom of the accumulator: 4096 rows of the largest double
    out.append([MAXF] * 2048 + [-MAXF] * 2047 + [TINY])
    return out


def expect_cross(pairs):
    s = sum((Fraction(x) * Fraction(y) for x, y in pairs), Fraction(0))
    try:
        return "%016x 0" % struct.unpack("<Q", struct.pack("<d", float(s)))[0]
    except OverflowError:
        return "%016x 1" % 0


def cross_families(r, count):
    """vectors of PAIRS (x, y): the families of measure_acc.py against each other (so exponents of both signs and every magnitude meet), and the extremes of the product"""
    gens = [ma.fam_uniform, ma.mk(-1074, 1023), ma.mk(400, 1023), ma.mk(-1074, -480), ma.mk(-300, 300), ma.mk(-60, 60), ma.fam_cancel, ma.fam_const, ma.mk(-1074, -1023), ma.fam_edge]
    out = []
    for i in range(count):
        n = [1, 2, 3, 7, 30, 120][i % 6]
        xs = gens[i % len(gens)](r, n)[0]
        ys = gens[(i * 7 + 3) % len(gens)](r, n)[0]
        m = min(len(xs), len(ys))
        out.append(list(zip(xs[:m], ys[:m])))
    out.append([(MAXF, MAXF)])
    out.append([(MAXF, 1.0), (-MAXF, 1.0), (MAXF, 0.5)])
    out.append([(TINY, TINY)] * 7)
    out.append([(TINY, 1.0), (-TINY, 1.0), (TINY, 3.0)])
    out.append([(2.0 ** 600, 2.0 ** 600), (2.0 ** 600, -(2.0 ** 600))])
    out.append([(1e150, 1e-150), (-1e150, 1e-150), (1e150, -1e-150)])
    out.append([(float((-1) ** i) * 3.0, 1.0 + i) for i in range(100)])
    out.append([(MAXF, MAXF)] * 3 + [(-MAXF, MAXF)] * 3 + [(TINY, TINY)])
    out.append([(2.0 ** -537, 2.0 ** -537)] * 50)
    return out


def run(driver, jobs):
    inp = "".join("%d %d %s\n" % (k, len(xs), " ".join(hx(x) for x in xs)) for k, xs in jobs)
    p = subprocess.run([driver], input=inp, capture_output=True, text=True, timeout=1800)
    if p.returncode != 0:
        raise RuntimeError("driver exit %d: %s" % (p.returncode, p.stderr[-200:]))
    out = [ln.strip() for ln in p.stdout.strip().split("\n")]
    if len(out) != len(jobs):
        raise RuntimeError("%d answers for %d jobs" % (len(out), len(jobs)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", default=str(ROOT / "build" / "psum-check"))
    ap.add_argument("--vectors", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = SM64(a.seed)
    vecs = families(r, a.vectors)
    jobs, want = [], []
    for xs in vecs:
        for k in (1, 2, 3, 4):
            jobs.append((k, xs))
            want.append(expect(xs, k))
    cross = cross_families(r, a.vectors)
    for pairs in cross:
        jobs.append((5, [v for pr in pairs for v in pr]))
        want.append(expect_cross(pairs))
    shuf = []
    for k, xs in jobs:
        if k == 5:
            prs = [(xs[2 * i], xs[2 * i + 1]) for i in range(len(xs) // 2)]
            for i in range(len(prs) - 1, 0, -1):
                j = r.next() % (i + 1)
                prs[i], prs[j] = prs[j], prs[i]
            shuf.append((k, [v for pr in prs for v in pr]))
            continue
        ys = list(xs)
        for i in range(len(ys) - 1, 0, -1):
            j = r.next() % (i + 1)
            ys[i], ys[j] = ys[j], ys[i]
        shuf.append((k, ys))
    try:
        got = run(a.driver, jobs)
        got2 = run(a.driver, shuf)
    except Exception as e:
        print("gate_psum: harness failure: %s" % e)
        return 2
    bad = [(j, w, g) for j, w, g in zip(jobs, want, got) if w != g]
    order = [(j, g, g2) for j, g, g2 in zip(jobs, got, got2) if g != g2]
    flagged = sum(1 for w in want if w.endswith(" 1"))
    print("gate_psum: %d sums (%d vectors x k = 1..4 and %d vectors of products, %d beyond the largest double and flagged): %d differ from Fraction, %d differ after shuffling the rows" % (len(jobs), len(vecs), len(cross), flagged, len(bad), len(order)))
    for j, w, g in bad[:4]:
        print("  k=%d n=%d want %s got %s" % (j[0], len(j[1]), w, g))
    return 1 if (bad or order) else 0


if __name__ == "__main__":
    sys.exit(main())
