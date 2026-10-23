#!/usr/bin/env python3
"""S2.1 gate C-A (docs/corr.md section 5): the exact cross sums sum(x_i * y_i) of lib/acc/xsum.cho against `fractions.Fraction`, rounded once, on pairs of vectors drawn from the families of
measure_acc.py (uniform, every finite exponent, huge, tiny, subnormal, cancellation, constant, edges of the range, each family against every other) and the extremes. Also: the rows run again in a
shuffled order (the pairs move together) must print the same bytes; the roles of x and y swapped must print the same bytes (the symmetry r(x, y) = r(y, x) rests on it); and sum(x * x) must equal the
power sum k = 2 of lib/acc/psum.cho (the two accumulators share one unit, 2^-2148, which the coefficient of corr relies on).

    python3 scripts/gate_xsum.py [--driver build/xsum-check] [--psum build/psum-check] [--vectors 300] [--seed 1]

Exit 0: every answer is the exactly rounded value (or code 1 exactly when the exact sum is beyond the largest double) and every shuffle, swap and psum comparison agrees. Exit 1: a difference. Exit 2: the harness failed.
"""
import argparse
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


def expect(pairs):
    s = sum((Fraction(x) * Fraction(y) for x, y in pairs), Fraction(0))
    try:
        return "%016x 0" % struct.unpack("<Q", struct.pack("<d", float(s)))[0]
    except OverflowError:
        return "%016x 1" % 0


def families(r, count):
    gens = [ma.fam_uniform, ma.mk(-1074, 1023), ma.mk(400, 1023), ma.mk(-1074, -480), ma.mk(-300, 300), ma.mk(-60, 60), ma.fam_cancel, ma.fam_const, ma.mk(-1074, -1023), ma.fam_edge]
    out = []
    for i in range(count):
        n = [1, 2, 3, 7, 30, 200][i % 6]
        gx = gens[i % len(gens)]
        gy = gens[(i * 7 + i // len(gens)) % len(gens)]  # every family meets every other over the run
        xs = gx(r, n)[0]
        ys = gy(r, n)[0]
        out.append(list(zip(xs, ys[:len(xs)])) if len(ys) >= len(xs) else list(zip(xs[:len(ys)], ys)))
    out.append([(MAXF, MAXF)])
    out.append([(MAXF, -MAXF), (MAXF, MAXF)])
    out.append([(TINY, TINY)] * 20)
    out.append([(TINY, MAXF)] * 20)
    out.append([(2.0 ** 1000, 2.0 ** 1000)] * 2)
    out.append([(float((-1) ** i) * 3.0, 5.0) for i in range(100)])
    out.append([(1e150, 1e-150), (1e-150, -1e150), (1e150, 1e150)])
    out.append([(MAXF, MAXF)] * 257)    # just above 2^2106 units: a limb count without headroom wraps to a finite double
    out.append([(MAXF, MAXF)] * 4096)   # the headroom of the accumulator
    out.append([(MAXF, MAXF)] * 2048 + [(MAXF, -MAXF)] * 2047 + [(TINY, TINY)])
    out.append([(0.0, 5.0), (-0.0, 7.0), (3.0, 0.0)])  # zeros of both signs contribute nothing
    return out


def run(driver, vecs, square=False):
    inp = "".join("%d %s\n" % (len(v), " ".join("%s %s" % (hx(x), hx(y)) for x, y in v)) for v in vecs)
    p = subprocess.run([driver], input=inp, capture_output=True, text=True, timeout=1800)
    if p.returncode != 0:
        raise RuntimeError("driver exit %d: %s" % (p.returncode, p.stderr[-200:]))
    out = [ln.strip() for ln in p.stdout.strip().split("\n")]
    if len(out) != len(vecs):
        raise RuntimeError("%d answers for %d vectors" % (len(out), len(vecs)))
    return out


def run_psum2(driver, vecs):
    inp = "".join("2 %d %s\n" % (len(v), " ".join(hx(x) for x, _ in v)) for v in vecs)
    p = subprocess.run([driver], input=inp, capture_output=True, text=True, timeout=1800)
    if p.returncode != 0:
        raise RuntimeError("psum driver exit %d: %s" % (p.returncode, p.stderr[-200:]))
    return [ln.strip() for ln in p.stdout.strip().split("\n")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", default=str(ROOT / "build" / "xsum-check"))
    ap.add_argument("--psum", default=str(ROOT / "build" / "psum-check"))
    ap.add_argument("--vectors", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = SM64(a.seed)
    vecs = families(r, a.vectors)
    want = [expect(v) for v in vecs]
    shuf, swap = [], []
    for v in vecs:
        w = list(v)
        for i in range(len(w) - 1, 0, -1):
            j = r.next() % (i + 1)
            w[i], w[j] = w[j], w[i]
        shuf.append(w)
        swap.append([(y, x) for x, y in v])
    squares = [[(x, x) for x, _ in v] for v in vecs]
    try:
        got = run(a.driver, vecs)
        got_shuf = run(a.driver, shuf)
        got_swap = run(a.driver, swap)
        got_sq = run(a.driver, squares)
        ps = run_psum2(a.psum, vecs)
    except Exception as e:
        print("gate_xsum: harness failure: %s" % e)
        return 2
    bad = [(v, w, g) for v, w, g in zip(vecs, want, got) if w != g]
    order = [v for v, g, g2 in zip(vecs, got, got_shuf) if g != g2]
    sym = [v for v, g, g2 in zip(vecs, got, got_swap) if g != g2]
    same = [v for v, g, p in zip(vecs, got_sq, ps) if g != p]
    flagged = sum(1 for w in want if w.endswith(" 1"))
    print("gate_xsum: %d sums (%d beyond the largest double and flagged): %d differ from Fraction, %d differ after shuffling the rows, %d differ with x and y swapped, %d differ from psum k = 2" % (len(vecs), flagged, len(bad), len(order), len(sym), len(same)))
    for v, w, g in bad[:4]:
        print("  n=%d want %s got %s" % (len(v), w, g))
    return 1 if (bad or order or sym or same) else 0


if __name__ == "__main__":
    sys.exit(main())
