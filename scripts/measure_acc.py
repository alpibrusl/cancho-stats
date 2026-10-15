#!/usr/bin/env python3
"""G4 / section 5 of docs/numerics.md: how wide must the accumulator of sum(x^2) and sum(x*y) be?

    python3 scripts/measure_acc.py [--bin build/acc-check] [--seed 1] [--vectors 300] [--report]

Per family of vectors it computes sum(x^2) and sum(x*y) EXACTLY (python `fractions.Fraction`) and rounds once (Fraction -> float is correctly rounded), and compares:
  A  two-product into cancho-table's `facc`, rounded once (acc.xprod + facc)
  B  the pinned wide accumulator `ml.fsq` (4196 bits), squares only
  C  the plain double sum in row order (what a streaming mean/variance update amounts to at its best)
Reported per family: vectors, and for A, B, C the number exact to the bit, and the largest relative error of the rest. Then the exactness BOUNDARY of A: for single products
x*x at every binary exponent of x, the fraction of 20 random mantissas for which A is exact, so the range rule R of the document is read off a measurement.
Exit 1 if rule R (|x| zero or in [2^-485, 2^501)) admits an inexact A on any vector inside it (the gate G4), or if B is inexact anywhere (squares of finite doubles).
"""
import argparse
import math
import pathlib
import struct
import subprocess
import sys
from fractions import Fraction

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from gate_special import SM64, hx, unhx  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
MAXF = struct.unpack("<d", struct.pack("<Q", 0x7FEFFFFFFFFFFFFF))[0]


def rounded(fr):
    try:
        return float(fr)
    except OverflowError:
        return math.inf if fr > 0 else -math.inf


def rand_double(r, elo, ehi, signs=True):
    """a finite double with binary exponent uniform in [elo, ehi] (the exponent of the leading bit; elo >= -1074)."""
    e = elo + r.next() % (ehi - elo + 1)
    if e >= -1022:
        m = (1 << 52) | (r.next() & ((1 << 52) - 1))
        bits = ((e + 1023) << 52) | (m & ((1 << 52) - 1))
    else:
        # subnormal: leading bit at position e + 1074
        k = e + 1074
        m = (1 << k) | (r.next() & ((1 << k) - 1)) if k > 0 else 1
        bits = m
    if signs and (r.next() & 1):
        bits |= 1 << 63
    return struct.unpack("<d", struct.pack("<Q", bits))[0]


def fam_uniform(r, n):
    return [r.u() for _ in range(n)], [r.u() for _ in range(n)]


def mk(elo, ehi, signs=True):
    def f(r, n):
        return [rand_double(r, elo, ehi, signs) for _ in range(n)], [rand_double(r, elo, ehi, signs) for _ in range(n)]
    return f


def fam_cancel(r, n):
    xs, ys = [], []
    for i in range(n):
        v = rand_double(r, -20, 20)
        k = i % 3
        xs.append(v if k != 1 else -v)
        ys.append(rand_double(r, -20, 20) if k != 2 else 1e-30)
    return xs, ys


def fam_const(r, n):
    v = rand_double(r, -10, 10, signs=False)
    return [v] * n, [v] * n


def fam_edge(r, n):
    return [r.pick([2.0 ** -485, 2.0 ** 500 * 1.9999999999999998, 2.0 ** -486, 2.0 ** 501, 5e-324, 2.0 ** -1022, MAXF / 2 ** 513]) for _ in range(n)], \
           [r.pick([2.0 ** -485, 2.0 ** 500, 1.0, 3.0]) for _ in range(n)]


FAMILIES = [
    ("uniform [0,1)", fam_uniform, True),
    ("wide exponent, all finite", mk(-1074, 1023), False),
    ("large 2^[400, 510]", mk(400, 510), False),
    ("small 2^[-1074, -480]", mk(-1074, -480), False),
    ("inside rule R, 2^[-485, 500]", mk(-485, 500), True),
    ("moderate 2^[-60, 60]", mk(-60, 60), True),
    ("cancellation (x, -x, 1e-30)", fam_cancel, True),
    ("constant column", fam_const, True),
    ("subnormals only", mk(-1074, -1023), False),
    ("edges of rule R", fam_edge, None),
]


def in_rule(v):
    a = abs(v)
    return a == 0.0 or (2.0 ** -485 <= a < 2.0 ** 501)


def run(binary, vecs):
    lines = []
    for xs, ys in vecs:
        lines.append("%d %s %s" % (len(xs), " ".join(hx(x) for x in xs), " ".join(hx(y) for y in ys)))
    r = subprocess.run([binary], input="\n".join(lines) + "\n", capture_output=True, text=True, timeout=1200)
    if r.returncode != 0:
        raise RuntimeError("acc-check exit %d: %s" % (r.returncode, r.stderr[-300:]))
    out = [ln.split() for ln in r.stdout.strip().split("\n")]
    if len(out) != len(vecs):
        raise RuntimeError("acc-check answered %d lines for %d vectors" % (len(out), len(vecs)))
    return out


def relerr(got, want):
    if got == want:
        return 0.0
    if want == 0 or math.isinf(want) or math.isinf(got) or got != got:
        return math.inf
    return abs(got - want) / abs(want)


def exact_sums(xs, ys):
    sq = sum((Fraction(x) * Fraction(x) for x in xs), Fraction(0))
    xy = sum((Fraction(x) * Fraction(y) for x, y in zip(xs, ys)), Fraction(0))
    return rounded(sq), rounded(xy)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "acc-check"))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--vectors", type=int, default=300)
    a = ap.parse_args()
    r = SM64(a.seed)
    fail = []
    print("%-30s %6s | %-24s | %-24s | %-24s" % ("family", "n vec", "A sq: exact / worst rel", "A xy: exact / worst rel", "B sq: exact / worst rel"))
    print("%-30s %6s | %-24s | %-24s" % ("", "", "C sq (plain double)", "C xy (plain double)"))
    for name, gen, inside in FAMILIES:
        vecs = []
        for i in range(a.vectors):
            n = [1, 2, 3, 10, 100, 1000][i % 6]
            vecs.append(gen(r, n))
        out = run(a.bin, vecs)
        ea = [0, 0, 0, 0, 0]  # exact counts A sq, A xy, B sq, C sq, C xy
        refused = [0, 0, 0]
        worst = [0.0] * 5
        inside_inexact = 0
        for (xs, ys), o in zip(vecs, out):
            wsq, wxy = exact_sums(xs, ys)
            got = [unhx(o[0]), unhx(o[1]), unhx(o[2]), unhx(o[3]), unhx(o[4])]
            want = [wsq, wxy, wsq, wsq, wxy]
            vec_in = all(in_rule(v) for v in xs) and all(in_rule(v) for v in ys)
            stat = [int(o[5]), int(o[6]), int(o[7])]
            for k in range(5):
                if k < 3 and math.isinf(want[k]) and stat[k] > 0:
                    refused[k] += 1  # the sum is beyond the largest double and the status says so: a refusal, not an error
                    continue
                e = relerr(got[k], want[k])
                if e == 0.0:
                    ea[k] += 1
                else:
                    worst[k] = max(worst[k], e)
                    if k in (0, 1) and vec_in:
                        inside_inexact += 1
                    if k == 2:
                        fail.append("B inexact in family %s" % name)
        print("%-30s %6d | %5d / %-16.3g | %5d / %-16.3g | %5d / %-16.3g   (status-refused A sq %d, A xy %d, B %d)" % (name, len(vecs), ea[0], worst[0], ea[1], worst[1], ea[2], worst[2], refused[0], refused[1], refused[2]))
        print("%-30s %6s | %5d / %-16.3g | %5d / %-16.3g" % ("", "", ea[3], worst[3], ea[4], worst[4]))
        if inside_inexact:
            fail.append("A inexact INSIDE rule R in family %s: %d results" % (name, inside_inexact))
    # The boundary of A, found by cancellation so that the low part of a product cannot hide under the final rounding: the pair (x, p) times (y, -1) has the exact sum
    # x*y - p, which is the rounding error of the product itself. A is exact there iff it kept the low part of x*y. (For s above 900 the pair is (p/256, -256), a power of two, so that p itself is never split.) s = ex + ey is the sum of the binary exponents.
    print("\nexactness boundary of A (sum x*y - fl(x*y), by cancellation), by s = ex + ey, 20 random mantissa pairs per s, three shapes d = ex - ey")
    window_bad = []
    for d in (-300, 0, 300):
        vecs, meta = [], []
        for s_ in list(range(-1100, -900)) + list(range(960, 1045)):
            ex = (s_ + d) // 2
            ey = s_ - ex
            if not (-1074 <= ex <= 1023 and -1074 <= ey <= 1023):
                continue
            for _ in range(20):
                x = rand_double(r, ex, ex, signs=False)
                y = rand_double(r, ey, ey, signs=False)
                p = x * y
                if math.isinf(p):
                    continue
                sc = 256.0 if s_ > 900 else 1.0  # p must not be split by the Dekker product itself (|p| < 2^996): an exact power-of-two scale
                vecs.append(([x, p / sc], [y, -sc]))
                meta.append(s_)
        out = run(a.bin, vecs)
        frac = {}
        for (xs, ys), o, s_ in zip(vecs, out, meta):
            want = rounded(Fraction(xs[0]) * Fraction(ys[0]) + Fraction(xs[1]) * Fraction(ys[1]))
            got = unhx(o[1])
            frac.setdefault(s_, []).append(got == want and int(o[6]) == 0)
        prev = None
        line = []
        for s_ in sorted(frac):
            k = sum(frac[s_])
            n_ = len(frac[s_])
            state = "all" if k == n_ else ("none" if k == 0 else "some")
            if state != prev:
                line.append("s=%d:%s" % (s_, state))
                prev = state
        print("  d=%4d  %s" % (d, "  ".join(line)))
        for s_ in range(-970, 1001):
            if s_ in frac and not all(frac[s_]):
                window_bad.append((d, s_))
    # the splitting overflow: |x| >= 2^996 with y = 1/x so that the product is moderate
    vecs, meta = [], []
    for ex in range(985, 1016):
        for _ in range(20):
            x = rand_double(r, ex, ex, signs=False)
            y = 1.0 / x
            p = x * y
            vecs.append(([x, p], [y, -1.0]))
            meta.append(ex)
    out = run(a.bin, vecs)
    frac = {}
    for (xs, ys), o, ex in zip(vecs, out, meta):
        want = rounded(Fraction(xs[0]) * Fraction(ys[0]) - Fraction(xs[1]))
        frac.setdefault(ex, []).append(unhx(o[1]) == want and int(o[6]) == 0)
    prev = None
    line = []
    for ex in sorted(frac):
        k = sum(frac[ex])
        state = "all" if k == 20 else ("none" if k == 0 else "some")
        if state != prev:
            line.append("ex=%d:%s" % (ex, state))
            prev = state
    print("  splitting overflow (x = 2^ex.., y = 1/x):  " + "  ".join(line))
    if window_bad:
        fail.append("A is inexact inside the window s in [-970, 1000]: %s" % window_bad[:4])
    print()
    if fail:
        print("measure_acc: FAILED: " + "; ".join(sorted(set(fail))))
        return 1
    print("measure_acc: rule R holds on every vector sampled; B exact everywhere")
    return 0


if __name__ == "__main__":
    sys.exit(main())
