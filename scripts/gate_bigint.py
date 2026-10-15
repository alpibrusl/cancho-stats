#!/usr/bin/env python3
"""S1-A part 1 (docs/describe.md section 5): the big-integer arithmetic of lib/acc/bigint.cho against Python's own integers and fractions.

    python3 scripts/gate_bigint.py [--driver build/bigint-check] [--n 3000] [--seed 1]

Cases: products, sums and differences of random numbers of 1 to 4,700 bits; ratios `(N / D) * 2^E` rounded to a double, against `float(Fraction(N, D) * 2^E)` (Python's conversion is correctly
rounded, round half to even, subnormals included); and the hard ratios: exact ties at 53 bits, one above and one below a tie, the subnormal boundaries, the largest double and the first value beyond
it (code 1), powers of two, and numbers with a 1 followed by 60 zeros. Exit 0: every answer equals Python's. Exit 1: a difference. Exit 2: the harness failed.
"""
import argparse
import pathlib
import struct
import subprocess
import sys
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402

MAXF = float.fromhex("0x1.fffffffffffffp+1023")


def rnd_int(r, bits):
    if bits <= 0:
        return 0
    v = 0
    for _ in range((bits + 63) // 64):
        v = (v << 64) | r.next()
    return v & ((1 << bits) - 1) | (1 << (bits - 1))


def bits_of(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def expect_ratio(n, d, e):
    v = Fraction(n, d) * (Fraction(2) ** e)
    try:
        return "%016x 0" % bits_of(float(v))
    except OverflowError:
        return "%016x 1" % 0


def cases(r, count):
    out = []
    for i in range(count):
        k = i % 5
        if k == 0:
            a, b = rnd_int(r, 1 + r.next() % 4700), rnd_int(r, 1 + r.next() % 4700)
            out.append(("mul %x %x" % (a, b), "%x" % (a * b)))
        elif k == 1 and i % 10 == 1:
            a, b = rnd_int(r, 1 + r.next() % 4700), rnd_int(r, 1 + r.next() % 4700)
            out.append(("add %x %x" % (a, b), "%x" % (a + b)))
        elif k == 1:
            a, b = rnd_int(r, 1 + r.next() % 4700), rnd_int(r, 1 + r.next() % 4700)
            if a < b:
                a, b = b, a
            out.append(("sub %x %x" % (a, b), "%x" % (a - b)))
        else:
            n = rnd_int(r, 1 + r.next() % 4700)
            d = rnd_int(r, 1 + r.next() % 130)
            e = int(r.next() % 2300) - 1150
            out.append(("ratio %x %x %d" % (n, d, e), expect_ratio(n, d, e)))
    # the hard ones
    hard = []
    for t in range(60):
        base = (1 << 53) + 1  # a tie after one bit
        hard.append((base << t, 1, 0))
        hard.append(((base << t) + 1, 1, 0))
        hard.append(((base << t) - 1, 1, 0))
        hard.append(((1 << 54) + 3, 1 << 1, t))
    for e in (-1080, -1075, -1074, -1073, -1023, -1022, -1021, 970, 971, 972, 1020, 1023, 1024):
        for n in (1, 3, (1 << 52) + 1, (1 << 53) - 1, (1 << 53) + 1, (1 << 53) + 2, (1 << 54) - 1):
            hard.append((n, 1, e))
    hard.append(((1 << 1024) - (1 << 970), 1, 0))  # the largest double
    hard.append(((1 << 1024) - (1 << 970) + (1 << 969), 1, 0))  # the tie above it: rounds to infinity
    hard.append(((1 << 1024) - (1 << 970) + (1 << 969) - 1, 1, 0))  # one below the tie
    for dd in (3, 7, 10 ** 20 + 7, (1 << 70) + 1):
        for tie in ((1 << 53) + 1, (1 << 53) + 3, (1 << 52) + 1):
            for t in (0, 5, -1074 + 1):
                hard.append((tie * dd, dd, t))
                hard.append((tie * dd + 1, dd, t))
                hard.append((tie * dd - 1, dd, t))
    hard.append((1, 1 << 60, 0))
    hard.append(((1 << 100) + 1, 1 << 100, 0))
    for n, d, e in hard:
        out.append(("ratio %x %x %d" % (n, d, e), expect_ratio(n, d, e)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", default=str(ROOT / "build" / "bigint-check"))
    ap.add_argument("--n", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = SM64(a.seed)
    cs = cases(r, a.n)
    inp = "\n".join(c[0] for c in cs) + "\n"
    try:
        p = subprocess.run([a.driver], input=inp, capture_output=True, text=True, timeout=1800)
    except subprocess.TimeoutExpired:
        print("gate_bigint: the driver did not finish (hang)")
        return 2
    if p.returncode != 0:
        print("gate_bigint: the driver exited %d (trap): %s" % (p.returncode, p.stderr[-200:]))
        return 2
    got = p.stdout.strip().split("\n")
    if len(got) != len(cs):
        print("gate_bigint: %d answers for %d cases" % (len(got), len(cs)))
        return 2
    bad = [(c, g) for c, g in zip(cs, got) if c[1] != g.strip()]
    kinds = {}
    for c in cs:
        kinds[c[0].split()[0]] = kinds.get(c[0].split()[0], 0) + 1
    print("gate_bigint: %d cases (%s): %d differ" % (len(cs), ", ".join("%s %d" % kv for kv in sorted(kinds.items())), len(bad)))
    for c, g in bad[:5]:
        print("  %s...\n    want %s\n    got  %s" % (c[0][:90], c[1][:60], g[:60]))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
