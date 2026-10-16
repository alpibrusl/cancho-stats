#!/usr/bin/env python3
"""S1-B (docs/describe.md section 5), for the statistics `stats describe` prints so far: count, missing, kind, min, max, mean, variance, sd, skewness and excess kurtosis of generated CSV files against `fractions.Fraction`.

    python3 scripts/gate_describe.py [--bin build/stats] [--files 300] [--seed 1]

Each file has columns of the shapes that break a float sum (a large offset with small noise, wide exponents, huge and tiny values, constants, one row, all cells missing, text, quoted numbers, signs), written as
decimal text. The expected mean and variance are the exact rational values of the doubles the text reads as, rounded ONCE by Python's Fraction -> float (correctly rounded); sd is the square root of that
rounded variance (math.sqrt is correctly rounded: the two roundings the answer states). The comparison is on the doubles parsed from the printed JSON, bit for bit. Every file is also run with its rows shuffled,
and with --ddof 0, and the shuffled answer must equal the original byte for byte. Variance is null exactly when n - ddof < 1 or the exact variance is beyond the largest double.
Exit 0: every number is the expected one. Exit 1: a difference. Exit 2: the harness failed.
"""
import argparse
import json
import math
import pathlib
import struct
import subprocess
import sys
import tempfile
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
import gate_quantile  # noqa: E402

MAXF = float.fromhex("0x1.fffffffffffffp+1023")


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def to_float(fr):
    try:
        return float(fr)
    except OverflowError:
        return None


def column(r, shape, n):
    """a list of cell strings (the text of the CSV cells) for one column"""
    cells = []
    for i in range(n):
        if shape == "uniform":
            cells.append(repr(r.u()))
        elif shape == "ints":
            cells.append(str(int(r.next() % 2001) - 1000))
        elif shape == "offset":
            cells.append(repr(1e8 + r.u()))
        elif shape == "offset-neg":
            cells.append(repr(-123456789.0 + r.u() * 1e-3))
        elif shape == "wide":
            e = int(r.next() % 600) - 300
            cells.append(repr(math.ldexp(1.0 + r.u(), e) * (1 if r.next() & 1 else -1)))
        elif shape == "tiny":
            cells.append(repr(math.ldexp(1.0 + r.u(), -1000 - int(r.next() % 60))))
        elif shape == "huge":
            cells.append(repr(math.ldexp(1.0 + r.u(), 600 + int(r.next() % 400)) * (1 if r.next() & 1 else -1)))
        elif shape == "const":
            cells.append("0.1")
        elif shape == "gappy":
            cells.append("" if r.next() % 3 == 0 else repr(r.u() * 100))
        elif shape == "empty":
            cells.append("")
        elif shape == "text":
            cells.append(r.pick(["Paris", "Rome", "Oslo", "x1"]))
        elif shape == "mixed":
            cells.append(r.pick(["1", "2.5", "n/a"]))
        elif shape == "quoted":
            cells.append('"%s"' % repr(r.u()))
        elif shape == "exp":
            cells.append("%de%d" % (int(r.next() % 9) + 1, int(r.next() % 40) - 20))
        else:
            raise KeyError(shape)
    return cells


SHAPES = ["uniform", "ints", "offset", "offset-neg", "wide", "tiny", "huge", "const", "gappy", "empty", "text", "mixed", "quoted", "exp"]


def sqrt_round(x):
    """the correctly rounded square root of the non-negative Fraction x: the integer root of x 4^k with a sticky bit below the guard bit"""
    if x == 0:
        return 0.0
    num, den = x.numerator, x.denominator
    k = max(0, (132 + den.bit_length() - num.bit_length() + 1) // 2)
    y = (num << (2 * k)) // den
    sroot = math.isqrt(y)
    inexact = 1 if sroot * sroot * den != (num << (2 * k)) else 0
    return float(Fraction(sroot * 2 + inexact, 2 ** (k + 1)))


def shape_moments(fr):
    """(skewness, excess kurtosis) of the exact values, each rounded once, or (None, None) when the variance is zero"""
    n = len(fr)
    mu = sum(fr, Fraction(0)) / n
    d = [x - mu for x in fr]
    m2 = sum((x * x for x in d), Fraction(0)) / n
    if m2 == 0:
        return None, None
    m3 = sum((x ** 3 for x in d), Fraction(0)) / n
    m4 = sum((x ** 4 for x in d), Fraction(0)) / n
    g1 = sqrt_round(m3 * m3 / (m2 ** 3))
    if m3 < 0:
        g1 = -g1
    return g1, float(m4 / (m2 * m2) - 3)


def expected(cells, ddof, qtype=7, probs=("0", "0.25", "0.5", "0.75", "1")):
    vals = []
    kind = "numeric"
    nonmiss = 0
    for c in cells:
        if c == "":
            continue
        nonmiss += 1
        if c.startswith('"'):
            kind = "text"
            continue
        try:
            v = float(c)
            if math.isinf(v) or math.isnan(v):
                raise ValueError
            vals.append(v)
        except ValueError:
            kind = "text"
    out = {"count": nonmiss, "missing": len(cells) - nonmiss}
    if nonmiss == 0:
        out["kind"] = "empty"
        return out
    out["kind"] = kind
    if kind == "text":
        return out
    n = len(vals)
    fr = [Fraction(v) for v in vals]
    s1 = sum(fr, Fraction(0))
    s2 = sum((x * x for x in fr), Fraction(0))
    out["min"] = min(vals)
    out["max"] = max(vals)
    out["mean"] = float(s1 / n)
    if n - ddof >= 1:
        var = (n * s2 - s1 * s1) / (n * (n - ddof))
        v = to_float(var)
        out["variance"] = v
        out["sd"] = math.sqrt(v) if v is not None else None
    else:
        out["variance"] = None
        out["sd"] = None
    srt = sorted(vals)
    out["quantiles"] = [float(gate_quantile.hf(srt, Fraction(pt), qtype)) for pt in probs]
    if n >= 2:
        out["skewness"], out["kurtosis"] = shape_moments(fr)
    else:
        out["skewness"], out["kurtosis"] = None, None
    return out


def run(binary, tmp, names, rows, ddof, extra=()):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    p = subprocess.run([binary, "describe", "--root", str(tmp), "--ddof", str(ddof)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=120)
    return p


def same(a, b):
    if a is None or b is None:
        return a is b
    if isinstance(a, float) or isinstance(b, float):
        return bits(float(a)) == bits(float(b))
    return a == b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    checked = 0
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            ncols = 1 + int(r.next() % 5)
            shapes = [SHAPES[int(r.next() % len(SHAPES))] for _ in range(ncols)]
            n = [1, 2, 3, 5, 20, 100, 400][int(r.next() % 7)]
            cols = [column(r, s, n) for s in shapes]
            names = ["c%d_%s" % (j, s) for j, s in enumerate(shapes)]
            rows = [[cols[j][i] for j in range(ncols)] for i in range(n)]
            for ddof in (1, 0):
                try:
                    p = run(a.bin, tmp, names, rows, ddof)
                except subprocess.TimeoutExpired:
                    print("gate_describe: the program did not finish (hang)")
                    return 2
                if p.returncode != 0:
                    print("gate_describe: exit %d on file %d: %s" % (p.returncode, fi, p.stdout[:200]))
                    return 2
                doc = json.loads(p.stdout)
                got = doc["data"]["columns"]
                if doc["data"]["rows"] != n:
                    bad.append((fi, "rows"))
                for j in range(ncols):
                    want = expected(cols[j], ddof)
                    for k, wv in want.items():
                        checked += 1
                        if k == "quantiles":
                            gq = [q["value"] for q in got[j].get("quantiles", [])]
                            if len(gq) != len(wv) or any(not same(x, y) for x, y in zip(gq, wv)):
                                bad.append((fi, "%s ddof %d quantiles: want %r got %r" % (names[j], ddof, wv, gq)))
                            continue
                        if not same(got[j].get(k), wv):
                            bad.append((fi, "%s ddof %d key %s: want %r got %r" % (names[j], ddof, k, wv, got[j].get(k))))
                    for k in ("min", "max", "mean", "variance", "sd", "skewness", "kurtosis", "quantiles"):
                        if k not in want and k in got[j]:
                            bad.append((fi, "%s: unexpected key %s" % (names[j], k)))
                if ddof == 1:
                    qt = 1 + int(r.next() % 9)
                    probs = []
                    for _ in range(1 + int(r.next() % 4)):
                        d = int(r.next() % 10)
                        probs.append("%.*f" % (d, (int(r.next() % (10 ** d + 1))) / (10 ** d)) if d else str(int(r.next() % 2)))
                    p3 = run(a.bin, tmp, names, rows, ddof, ["--quantile-type", str(qt), "--quantiles", ",".join(probs)])
                    if p3.returncode != 0:
                        bad.append((fi, "quantile flags: exit %d: %s" % (p3.returncode, p3.stdout[:120])))
                    else:
                        got3 = json.loads(p3.stdout)["data"]["columns"]
                        for j in range(ncols):
                            want3 = expected(cols[j], ddof, qt, tuple(probs)).get("quantiles")
                            if want3 is not None:
                                checked += 1
                                gq = [q["value"] for q in got3[j].get("quantiles", [])]
                                if len(gq) != len(want3) or any(not same(x, y) for x, y in zip(gq, want3)):
                                    bad.append((fi, "%s type %d p %s: want %r got %r" % (names[j], qt, probs, want3, gq)))
                    # row order cannot change a byte
                    order = list(range(n))
                    for i in range(n - 1, 0, -1):
                        k = int(r.next() % (i + 1))
                        order[i], order[k] = order[k], order[i]
                    p2 = run(a.bin, tmp, names, [rows[i] for i in order], ddof)
                    if p2.stdout != p.stdout:
                        bad.append((fi, "shuffled rows give other bytes"))
    print("gate_describe: %d files, %d values compared, %d differ" % (a.files, checked, len(bad)))
    for b in bad[:6]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
