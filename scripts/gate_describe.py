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
        elif shape == "fifth":
            # one zero in five, the rest ones: m2 = (2/5)^2 and skewness exactly -1.5 (the root of a perfect square)
            cells.append("0" if i % 5 == 0 else "1")
        elif shape == "symm":
            # values in pairs v, -v: the third central moment is exactly 0 and the variance is not
            v = repr(float(int(r.next() % 50) + 1))
            cells.append(v if i % 2 == 0 else "-" + cells[i - 1])
        elif shape == "levels":
            # a skewed set of about 30 levels, some missing: ties and a long tail
            k = int(r.next() % 30)
            cells.append("" if r.next() % 7 == 0 else "k%d" % (k * k % 29 if r.next() & 1 else k % 4))
        elif shape == "qtext":
            cells.append(r.pick(['"a b"', '"say ""hi"""', '"x,y"', '""', "plain", '"plain"']))
        elif shape == "utf8":
            cells.append(r.pick(["caf\u00e9", "\u00fcber", "\u65e5\u672c", "Z", "z", "a"]))
        elif shape == "uniq":
            cells.append("id%07d" % (int(r.next() % 10 ** 7)))
        elif shape == "exp":
            cells.append("%de%d" % (int(r.next() % 9) + 1, int(r.next() % 40) - 20))
        else:
            raise KeyError(shape)
    return cells


SHAPES = ["uniform", "ints", "offset", "offset-neg", "wide", "tiny", "huge", "const", "gappy", "empty", "text", "mixed", "quoted", "exp", "levels", "qtext", "utf8", "uniq", "fifth", "symm"]


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


def histogram(vals, bins):
    """(bins, edges, counts) with the edges rounded once and the counts by exact comparison against the rational edges, the last bin closed"""
    fr = [Fraction(v) for v in vals]
    lo, hi = min(fr), max(fr)
    if lo == hi:
        return 1, [float(lo), float(hi)], [len(vals)]
    edges = [lo + (hi - lo) * i / bins for i in range(bins + 1)]
    counts = [0] * bins
    for x in fr:
        idx = sum(1 for i in range(1, bins) if x >= edges[i])
        counts[idx] += 1
    return bins, [float(e) for e in edges], counts


def expected(cells, ddof, qtype=7, probs=("0", "0.25", "0.5", "0.75", "1"), bins=10):
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
        counts = {}
        for c in cells:
            if c == "":
                continue
            lv = c[1:-1] if c.startswith('"') else c
            counts[lv] = counts.get(lv, 0) + 1
        out["distinct"] = len(counts)
        out["levels_all"] = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0].encode("utf-8")))
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
    out["distinct"] = len(set(vals))
    if bins > 0:
        out["histogram"] = histogram(vals, bins)
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
                        if k == "histogram":
                            gh = got[j].get("histogram", {})
                            if gh.get("bins") != wv[0] or gh.get("counts") != wv[2] or len(gh.get("edges", [])) != len(wv[1]) or any(not same(x, y) for x, y in zip(gh.get("edges", []), wv[1])):
                                bad.append((fi, "%s ddof %d histogram: want %r got %r" % (names[j], ddof, wv, gh)))
                            continue
                        if k == "quantiles":
                            gq = [q["value"] for q in got[j].get("quantiles", [])]
                            if len(gq) != len(wv) or any(not same(x, y) for x, y in zip(gq, wv)):
                                bad.append((fi, "%s ddof %d quantiles: want %r got %r" % (names[j], ddof, wv, gq)))
                            continue
                        if k == "levels_all":
                            gl = [(q["value"], q["count"]) for q in got[j].get("levels", [])]
                            if gl != wv[:10] or got[j].get("levels_capped") is not False:
                                bad.append((fi, "%s ddof %d levels: want %r got %r" % (names[j], ddof, wv[:10], gl)))
                            continue
                        if not same(got[j].get(k), wv):
                            bad.append((fi, "%s ddof %d key %s: want %r got %r" % (names[j], ddof, k, wv, got[j].get(k))))
                    for k in ("min", "max", "mean", "variance", "sd", "skewness", "kurtosis", "quantiles", "distinct", "histogram", "levels", "levels_capped"):
                        if k not in want and k in got[j] and not (k in ("levels", "levels_capped") and "levels_all" in want):
                            bad.append((fi, "%s: unexpected key %s" % (names[j], k)))
                if ddof == 1:
                    qt = 1 + int(r.next() % 9)
                    probs = []
                    for _ in range(1 + int(r.next() % 4)):
                        d = int(r.next() % 10)
                        probs.append("%.*f" % (d, (int(r.next() % (10 ** d + 1))) / (10 ** d)) if d else str(int(r.next() % 2)))
                    nbins = [0, 1, 2, 3, 5, 7, 20, 64][int(r.next() % 8)]
                    p3 = run(a.bin, tmp, names, rows, ddof, ["--quantile-type", str(qt), "--quantiles", ",".join(probs), "--bins", str(nbins)])
                    if p3.returncode != 0:
                        bad.append((fi, "quantile flags: exit %d: %s" % (p3.returncode, p3.stdout[:120])))
                    else:
                        got3 = json.loads(p3.stdout)["data"]["columns"]
                        for j in range(ncols):
                            w3 = expected(cols[j], ddof, qt, tuple(probs), nbins)
                            want3 = w3.get("quantiles")
                            if "histogram" in w3 or nbins == 0:
                                gh = got3[j].get("histogram")
                                wh = w3.get("histogram")
                                checked += 1
                                if wh is None:
                                    if gh is not None:
                                        bad.append((fi, "%s: histogram with --bins 0" % names[j]))
                                elif gh is None or gh.get("bins") != wh[0] or gh.get("counts") != wh[2] or any(not same(x, y) for x, y in zip(gh.get("edges", []), wh[1])):
                                    bad.append((fi, "%s bins %d: want %r got %r" % (names[j], nbins, wh, gh)))
                            if want3 is not None:
                                checked += 1
                                gq = [q["value"] for q in got3[j].get("quantiles", [])]
                                if len(gq) != len(want3) or any(not same(x, y) for x, y in zip(gq, want3)):
                                    bad.append((fi, "%s type %d p %s: want %r got %r" % (names[j], qt, probs, want3, gq)))
                    # --levels N: the first N levels, or none and no `distinct` for a text column with 0
                    nl = [0, 1, 2, 3, 10, 100][int(r.next() % 6)]
                    p4 = run(a.bin, tmp, names, rows, ddof, ["--levels", str(nl)])
                    if p4.returncode != 0:
                        bad.append((fi, "--levels %d: exit %d: %s" % (nl, p4.returncode, p4.stdout[:120])))
                    else:
                        got4 = json.loads(p4.stdout)["data"]["columns"]
                        for j in range(ncols):
                            w4 = expected(cols[j], ddof)
                            if w4["kind"] != "text":
                                continue
                            checked += 1
                            if nl == 0:
                                if "levels" in got4[j] or "distinct" in got4[j]:
                                    bad.append((fi, "%s: levels with --levels 0" % names[j]))
                            else:
                                gl = [(q["value"], q["count"]) for q in got4[j].get("levels", [])]
                                if gl != w4["levels_all"][:nl] or got4[j].get("distinct") != w4["distinct"]:
                                    bad.append((fi, "%s --levels %d: want %r got %r" % (names[j], nl, w4["levels_all"][:nl], gl)))
                    # --max-levels M: the pairs are added row by row, column by column; the column that would need one more is capped
                    ml = [1, 2, 3, 8, 30, 1000][int(r.next() % 6)]
                    p5 = run(a.bin, tmp, names, rows, ddof, ["--max-levels", str(ml)])
                    if p5.returncode != 0:
                        bad.append((fi, "--max-levels %d: exit %d: %s" % (ml, p5.returncode, p5.stdout[:120])))
                    else:
                        got5 = json.loads(p5.stdout)["data"]["columns"]
                        wants = [expected(cols[j], ddof) for j in range(ncols)]
                        seen, capped = set(), set()
                        for i in range(n):
                            for j in range(ncols):
                                c = cols[j][i]
                                if wants[j]["kind"] != "text" or c == "":
                                    continue
                                key = (j, c[1:-1] if c.startswith('"') else c)
                                if key not in seen:
                                    if len(seen) < ml:
                                        seen.add(key)
                                    else:
                                        capped.add(j)
                        for j in range(ncols):
                            if wants[j]["kind"] != "text":
                                continue
                            checked += 1
                            g = got5[j]
                            if j in capped:
                                if g.get("levels_capped") is not True or g.get("distinct") is not None or g.get("levels") is not None:
                                    bad.append((fi, "%s --max-levels %d: want capped, got %r" % (names[j], ml, g)))
                            else:
                                gl = [(q["value"], q["count"]) for q in g.get("levels", [])]
                                if g.get("levels_capped") is not False or g.get("distinct") != wants[j]["distinct"] or gl != wants[j]["levels_all"][:10]:
                                    bad.append((fi, "%s --max-levels %d: want exact levels, got %r" % (names[j], ml, g)))
                    # row order cannot change a byte
                    order = list(range(n))
                    for i in range(n - 1, 0, -1):
                        k = int(r.next() % (i + 1))
                        order[i], order[k] = order[k], order[i]
                    p2 = run(a.bin, tmp, names, [rows[i] for i in order], ddof)
                    if p2.stdout != p.stdout:
                        bad.append((fi, "shuffled rows give other bytes"))
    # fixed columns that reach the exact square root and the zero third moment on every run: one zero in five (skewness exactly -1.5), pairs v and -v (skewness 0)
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        fixed = {"fifth5": ["0", "1", "1", "1", "1"], "fifth20": ["0" if i % 5 == 0 else "1" for i in range(20)], "symm2": ["3", "-3"], "symm10": ["4", "-4", "7", "-7", "1", "-1", "2", "-2", "9", "-9"], "tenth": ["0" if i % 10 == 0 else "1" for i in range(10)]}
        names = list(fixed)
        rows = [[fixed[nm][i] if i < len(fixed[nm]) else "" for nm in names] for i in range(20)]
        p = run(a.bin, tmp, names, rows, 1)
        if p.returncode != 0:
            bad.append((-1, "fixed columns: exit %d" % p.returncode))
        else:
            got = json.loads(p.stdout)["data"]["columns"]
            for j, nm in enumerate(names):
                want = expected(fixed[nm], 1)
                for k in ("skewness", "kurtosis", "variance", "mean"):
                    checked += 1
                    if not same(got[j].get(k), want.get(k)):
                        bad.append((-1, "fixed %s %s: want %r got %r" % (nm, k, want.get(k), got[j].get(k))))
    # a file of one text column with blank lines: each is one missing cell, in both passes
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        (tmp / "f.csv").write_text("c\nx\n\ny\nx\n\n")
        p = subprocess.run([a.bin, "describe", "--root", str(tmp), "f.csv"], capture_output=True, text=True)
        col = json.loads(p.stdout)["data"]["columns"][0] if p.returncode == 0 else {}
        checked += 1
        if col.get("count") != 3 or col.get("missing") != 2 or [(q["value"], q["count"]) for q in col.get("levels", [])] != [("x", 2), ("y", 1)]:
            bad.append((-1, "one text column with blank lines: %r" % p.stdout[:200]))
    # more than 256 text columns: the column number of a key is two bytes, and the levels of column j and of column j + 256 must not meet
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        wide = 300
        names = ["w%d" % j for j in range(wide)]
        rows = [["t%d" % ((j * 7 + i) % 5) for j in range(wide)] for i in range(6)]
        p = run(a.bin, tmp, names, rows, 1, ["--max-cols", "400"])
        if p.returncode != 0:
            bad.append((-1, "wide file: exit %d: %s" % (p.returncode, p.stdout[:120])))
        else:
            got = json.loads(p.stdout)["data"]["columns"]
            for j in range(wide):
                want = expected([rows[i][j] for i in range(6)], 1)
                checked += 1
                gl = [(q["value"], q["count"]) for q in got[j].get("levels", [])]
                if gl != want["levels_all"][:10] or got[j].get("distinct") != want["distinct"]:
                    bad.append((-1, "wide column %d: want %r got %r" % (j, want["levels_all"][:10], gl)))
    print("gate_describe: %d files, %d values compared, %d differ" % (a.files, checked, len(bad)))
    for b in bad[:6]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
