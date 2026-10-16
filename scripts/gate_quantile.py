#!/usr/bin/env python3
"""S1-C (docs/describe.md section 5): the nine quantile definitions of lib/acc/quantile.cho against an independent implementation of Hyndman and Fan (1996), in exact rational arithmetic.

    python3 scripts/gate_quantile.py [--driver build/quantile-check] [--cases 6000] [--seed 1]

The reference is written here from the paper's formulas, not from any library: Q(p) = (1 - g) x[j] + g x[j + 1] with j = floor(n p + m) and g = n p + m - j; for types 1 to 3 the weight g is replaced by
the rule of the paper (type 1: 0 if g = 0 else 1; type 2: 1/2 if g = 0 else 1; type 3, m = -1/2: 0 if g = 0 and j is even else 1), for types 4 to 9 m is 0, 1/2, p, 1 - p, (p + 1)/3, p/4 + 3/8;
x[0] = x[1] and x[n + 1] = x[n]. The value is an exact Fraction rounded once (Fraction -> float is correctly rounded). The cancho side uses the (a, b) form nppm = a + p (n + 1 - a - b), so the two
share no algebra. Cases: random data of the families of measure_acc.py, ties and duplicates, one to 400 values, all nine types, p in {0, 1} and random decimals of up to 9 digits and the neighbours of the
jumps (p = j/n, (j - 1/2)/n, ...). Also reported, NOT gated: how many answers equal `numpy.quantile(method=...)` bit for bit (numpy computes in floats); no sentence says they "match".
Exit 0: every answer equals the reference bit for bit. Exit 1: a difference. Exit 2: the harness failed.
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


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def hf(xs_sorted, p, t):
    """the exact value (a Fraction) of the type-t quantile of the sorted floats at the Fraction p"""
    n = len(xs_sorted)

    def x(i):
        return Fraction(xs_sorted[min(max(i, 1), n) - 1])

    m = {1: Fraction(0), 2: Fraction(0), 3: Fraction(-1, 2), 4: Fraction(0), 5: Fraction(1, 2), 6: p, 7: 1 - p, 8: (p + 1) / 3, 9: p / 4 + Fraction(3, 8)}[t]
    v = n * p + m
    j = math.floor(v)
    g = v - j
    if t == 1:
        w = Fraction(0) if g == 0 else Fraction(1)
    elif t == 2:
        w = Fraction(1, 2) if g == 0 else Fraction(1)
    elif t == 3:
        w = Fraction(0) if (g == 0 and j % 2 == 0) else Fraction(1)
    else:
        w = g
    return (1 - w) * x(j) + w * x(j + 1)


NUMPY = {1: "inverted_cdf", 2: "averaged_inverted_cdf", 3: "closest_observation", 4: "interpolated_inverted_cdf", 5: "hazen", 6: "weibull", 7: "linear", 8: "median_unbiased", 9: "normal_unbiased"}


def random_p(r, n):
    k = r.next() % 5
    if k == 0:
        e = int(r.next() % 2)
        return Fraction(e), 0, e
    if k == 1:
        d = 1 + int(r.next() % 9)
        num = int(r.next() % (10 ** d + 1))
        return Fraction(num, 10 ** d), d, num
    if k == 2:
        # the jumps: p = j / n, (j - 1/2) / n, (j + 1/2) / n, rounded to 9 digits down and up
        j = int(r.next() % (n + 1))
        base = Fraction(j, n) if r.next() % 2 else Fraction(2 * j - 1, 2 * n) if j > 0 else Fraction(0)
        num = max(0, min(10 ** 9, math.floor(base * 10 ** 9) + int(r.next() % 3) - 1))
        return Fraction(num, 10 ** 9), 9, num
    if k == 3:
        d = 1 + int(r.next() % 3)
        num = int(r.next() % (10 ** d + 1))
        return Fraction(num, 10 ** d), d, num
    d = 2
    num = [25, 50, 75, 10, 90, 5, 95, 1, 99][int(r.next() % 9)]
    return Fraction(num, 100), d, num


def data(r, i):
    gens = [ma.fam_uniform, ma.mk(-30, 30), ma.mk(-1074, 1023), ma.fam_cancel, ma.fam_const, ma.mk(-300, 300), ma.mk(-1074, -1023)]
    g = gens[i % len(gens)]
    n = [1, 2, 3, 4, 5, 7, 10, 25, 100, 400][int(r.next() % 10)]
    xs = g(r, n)[0]
    if i % 5 == 0 and n > 2:  # ties
        xs = [xs[int(r.next() % max(1, n // 3))] for _ in range(n)]
    return xs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", default=str(ROOT / "build" / "quantile-check"))
    ap.add_argument("--cases", type=int, default=6000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = SM64(a.seed)
    jobs = []
    for i in range(a.cases):
        xs = data(r, i)
        t = 1 + i % 9
        p, k, num = random_p(r, len(xs))
        jobs.append((xs, t, p, k, num))
    inp = "".join("%d %d %d %d %s\n" % (t, num, k, len(xs), " ".join(hx(x) for x in xs)) for xs, t, p, k, num in jobs)
    try:
        pr = subprocess.run([a.driver], input=inp, capture_output=True, text=True, timeout=1800)
    except subprocess.TimeoutExpired:
        print("gate_quantile: the driver did not finish (hang)")
        return 2
    if pr.returncode != 0:
        print("gate_quantile: the driver exited %d (trap): %s" % (pr.returncode, pr.stderr[-200:]))
        return 2
    out = pr.stdout.strip().split("\n")
    if len(out) != len(jobs):
        print("gate_quantile: %d answers for %d cases" % (len(out), len(jobs)))
        return 2
    bad = []
    for (xs, t, p, k, num), line in zip(jobs, out):
        want = float(hf(sorted(xs), p, t))
        if want == 0.0:
            want = 0.0
        got = line.split()[0]
        if "%016x" % bits(want) != got:
            # zero has one spelling here: -0.0 and 0.0 are the same key
            if want == 0.0 and int(got, 16) in (0, 1 << 63):
                continue
            bad.append((t, len(xs), p, got, "%016x" % bits(want)))
    try:
        import numpy as np

        eq = 0
        tot = 0
        for (xs, t, p, k, num), line in zip(jobs[:1500], out):
            q = float(np.quantile(np.array(xs), float(p), method=NUMPY[t]))
            tot += 1
            if bits(q) == int(line.split()[0], 16) or (q == 0.0 and int(line.split()[0], 16) in (0, 1 << 63)):
                eq += 1
        note = "; reported only: %d of %d answers equal numpy.quantile(method=...) bit for bit (numpy rounds in floats)" % (eq, tot)
    except ImportError:
        note = ""
    print("gate_quantile: %d cases over all nine types: %d differ from the exact Hyndman-Fan reference%s" % (len(jobs), len(bad), note))
    for b in bad[:5]:
        print("  type %d n %d p %s: got %s want %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
