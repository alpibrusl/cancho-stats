#!/usr/bin/env python3
"""S2.2 gates C-B and C-C and C-D (docs/corr.md section 5), for `stats corr` (Pearson): generated CSV files against

    python3 scripts/gate_corr.py [--bin build/stats] [--files 200] [--seed 1] [--scipy]

Oracle, from `fractions.Fraction`, of the doubles the text reads as: for every pair of numeric columns (both modes `--missing pairwise` and `listwise`) the rows used, n, the exact integers
N = n Sxy - Sx Sy, Lx = n Sxx - Sx^2, Ly = n Syy - Sy^2; r = sign(N) sqrt(N^2 / (Lx Ly)) rounded ONCE (bytes compared); |r| = 1 decided by N^2 = Lx Ly; the statistic
t = sign(N) sqrt(N^2 (n - 2) / (Lx Ly - N^2)) rounded once (bytes compared), df = n - 2; and the two-sided p-value 2 tsf(|t|, n - 2) = I_{1 - r^2}((n - 2)/2, 1/2) from `mpmath` at 60 digits on the EXACT
rational 1 - r^2. The p-value bound is fixed here, before the first run (docs/corr.md section 5, C-C): relative 1e-12 + (n - 1) 2^-53 for p-values from 1e-300 up; for a smaller p the program must
answer less than 1e-290. A p-value is exactly 0 when |r| = 1. Undefined values are null with a reason: a constant column (Lx or Ly zero, which includes n < 2), n < 3 for the statistic and the p-value.
Also: the rows shuffled give the same bytes; `--cols` in another order gives the same bytes (the answer keeps the order of the file); the columns of the file reversed give the mirrored pairs with the same values; the three CONTROLS (n instead of n - 2 degrees of freedom, a one-sided tail, r computed
by a float Welford update) must fall outside the bound on at least a tenth of the cases, so the gate is able to fail. With --scipy: pearsonr agrees to 1e-12 absolute on r and 1e-9 relative on p on the
well-conditioned families (uniform, ints, lin) with |r| <= 0.999; the other pairs are reported, not gated (the first run showed scipy's p-value, which is computed from the rounded r with 1 - r^2 in floats, off by more than 1e-9 relative near |r| = 1, where this program agrees with the exact oracle).

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

import mpmath

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
import gate_describe as gd  # noqa: E402

mpmath.mp.dps = 60
EPS = 2.0 ** -53
GOOD = ("uniform", "ints", "lin")
SHAPES = ["uniform", "ints", "lin", "offset", "offset-neg", "wide", "tiny", "huge", "const", "gappy", "empty", "text", "mixed", "lin", "lin", "ints", "uniform"]


def _t(*cols):
    return [[repr(float(v)) for v in c] for c in cols]


# hand cases that every run includes: |t| beyond the largest double (n = 3, the tiny third row decides the last digits), r exactly -1 and +1, n = 2, subnormals, a huge exact pair
FIXED = [
    _t([1.0, 2.0, 2.0 ** -1000], [2.0, 4.0, 2 * 2.0 ** -1000 * (1 + 2.0 ** -52)]),
    _t([1, 2, 3, 4], [-2, -4, -6, -8]),
    _t([1, 2], [3, 5]),
    _t([1e300, 2e300, -3e300], [1e300, 2e300, -3e300]),
    _t([5e-324, 1e-323, 1.5e-323, 2e-323], [1e-323, 3e-323, 2e-323, 4e-323]),
    _t([1.0, 2.0, 3.0], [1.0, 2.0, 3.0000000000000004]),
    _t([1e8, 1e8 + 1, 1e8 + 2, 1e8 + 4], [3.0, 1.0, 2.0, 5.0], [7.0, 7.0, 7.0, 7.0]),
    _t([1, 2, 3, 4, 5], [1, 2, 3, 2, 1], [1, 2, 3, 2, 1.0000000000000002]),  # r = 0 exactly, and r about 1e-17: the p-value is 1, never above it
]


def column(r, shape, n, prev):
    """the cells of one column; `lin` follows the previous column (a correlated pair, with noise)"""
    if shape == "lin":
        out = []
        for i in range(n):
            try:
                base = float(prev[i])
            except (ValueError, IndexError):
                base = r.u()
            out.append(repr(base * 2.5 + (r.u() - 0.5) * 0.5))
        return out
    return gd.column(r, shape, n)


def parse(cell):
    """the double a cell reads as, or None for a missing cell; a ValueError for a cell that is not a number"""
    if cell == "":
        return None
    if cell.startswith('"'):
        raise ValueError("quoted")
    x = float(cell)
    if math.isnan(x) or math.isinf(x):
        raise ValueError("not finite")
    return x


def numeric_columns(cols):
    """the indices of the columns that are numeric: at least one number and every non-empty cell a number"""
    out = []
    for j, c in enumerate(cols):
        try:
            vals = [parse(x) for x in c]
        except ValueError:
            continue
        if any(v is not None for v in vals):
            out.append(j)
    return out


def exact_pair(xs, ys):
    """(n, N, Lx, Ly) of the paired values, exact integers/Fractions"""
    n = len(xs)
    fx = [Fraction(x) for x in xs]
    fy = [Fraction(y) for y in ys]
    sx = sum(fx, Fraction(0))
    sy = sum(fy, Fraction(0))
    sxx = sum((x * x for x in fx), Fraction(0))
    syy = sum((y * y for y in fy), Fraction(0))
    sxy = sum((x * y for x, y in zip(fx, fy)), Fraction(0))
    return n, n * sxy - sx * sy, n * sxx - sx * sx, n * syy - sy * sy


def expected(xs, ys):
    """the oracle of one pair: a dict with n, r, statistic, df, p (mpf or None), exact_one, reason"""
    n, N, Lx, Ly = exact_pair(xs, ys)
    out = {"n": n, "r": None, "statistic": None, "df": None, "p": None, "one": False, "reason": None}
    if Lx == 0 or Ly == 0:
        out["reason"] = "constant-column"
        return out
    D = Lx * Ly
    nn = N * N
    r = gd.sqrt_round(Fraction(nn) / D)
    if N < 0:
        r = -r
    out["r"] = r
    if n < 3:
        out["reason"] = "too-few-rows"
        return out
    out["df"] = n - 2
    if nn == D:
        out["one"] = True
        out["p"] = mpmath.mpf(0)
        return out
    try:
        t = gd.sqrt_round(Fraction(nn) * (n - 2) / (D - nn))
    except OverflowError:
        # |t| beyond the largest double: no statistic and no p-value (below 3.6e-309 for every df)
        out["reason"] = "statistic-beyond-double"
        return out
    out["statistic"] = -t if N < 0 else t
    x = Fraction(D - nn) / D
    out["p"] = mpmath.betainc(mpmath.mpf(n - 2) / 2, mpmath.mpf(1) / 2, 0, mpmath.mpf(x.numerator) / x.denominator, regularized=True)
    return out


def welford_r(xs, ys):
    """control: the float one-pass (Welford) coefficient"""
    n = 0
    mx = my = sxx = syy = sxy = 0.0
    for x, y in zip(xs, ys):
        n += 1
        dx = x - mx
        mx += dx / n
        dy = y - my
        my += dy / n
        sxx += dx * (x - mx)
        syy += dy * (y - my)
        sxy += dx * (y - my)
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def p_ok(got, want, n):
    """the bound of C-C"""
    if want == 0:
        return got == 0
    if want < mpmath.mpf("1e-300"):
        return got is not None and got < 1e-290
    return got is not None and abs(mpmath.mpf(got) - want) <= (mpmath.mpf("1e-12") + (n - 1) * EPS) * want


def run(binary, tmp, names, rows, extra=()):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "corr", "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scipy", action="store_true")
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    checked = 0
    ctl = {"df-n": [0, 0], "one-sided": [0, 0], "welford": [0, 0]}  # [outside the bound, cases]
    sp = {"r": 0.0, "p": 0.0, "n": 0, "skipped": 0}
    pst = {"max": 0.0, "ratio": 0.0, "n": 0}  # the largest relative error of a p-value against mpmath, and its largest fraction of the bound
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            if fi < len(FIXED):
                cols = FIXED[fi]
                shapes = ["fixed"] * len(cols)
                ncols = len(cols)
                n = len(cols[0])
            else:
                ncols = 2 + int(r.next() % 4)
                shapes = [SHAPES[int(r.next() % len(SHAPES))] for _ in range(ncols)]
                n = [1, 2, 3, 4, 5, 8, 20, 60, 200][int(r.next() % 9)]
                cols = []
                for j, s in enumerate(shapes):
                    cols.append(column(r, s, n, cols[j - 1] if j else []))
            names = ["c%d_%s" % (j, s) for j, s in enumerate(shapes)]
            rows = [[cols[j][i] for j in range(ncols)] for i in range(n)]
            num = numeric_columns(cols)
            for mode in ("pairwise", "listwise"):
                try:
                    p = run(a.bin, tmp, names, rows, ["--missing", mode])
                except subprocess.TimeoutExpired:
                    print("gate_corr: the program did not finish (hang)")
                    return 2
                if len(num) < 2:
                    ok = p.returncode != 0 and json.loads(p.stdout or p.stderr or "{}").get("error", {}).get("rule") == "stats.corr-columns"
                    checked += 1
                    if not ok:
                        bad.append((fi, "fewer than two numeric columns: want stats.corr-columns, got exit %d %s" % (p.returncode, p.stdout[:120])))
                    continue
                if p.returncode != 0:
                    keep = pathlib.Path(tempfile.gettempdir()) / ("gate_corr_failed_%d.csv" % fi)
                    keep.write_text((tmp / "f.csv").read_text())
                    print("gate_corr: exit %d on file %d (%s): %s (the file is %s)" % (p.returncode, fi, mode, p.stdout[:200], keep))
                    return 2
                doc = json.loads(p.stdout)["data"]
                vals = {j: [parse(x) for x in cols[j]] for j in num}
                if mode == "listwise":
                    used = [i for i in range(n) if all(vals[j][i] is not None for j in num)]
                want_pairs = [(x, y) for xi, x in enumerate(num) for y in num[xi + 1:]]
                got_pairs = doc["pairs"]
                if [(g["x"], g["y"]) for g in got_pairs] != [(names[x], names[y]) for x, y in want_pairs]:
                    bad.append((fi, "%s: the pairs are %r" % (mode, [(g["x"], g["y"]) for g in got_pairs])))
                    continue
                for g, (x, y) in zip(got_pairs, want_pairs):
                    if mode == "pairwise":
                        used = [i for i in range(n) if vals[x][i] is not None and vals[y][i] is not None]
                    xs = [vals[x][i] for i in used]
                    ys = [vals[y][i] for i in used]
                    w = expected(xs, ys)
                    tag = "%s %s~%s" % (mode, names[x], names[y])
                    checked += 1
                    if g["n"] != w["n"]:
                        bad.append((fi, "%s n: want %r got %r" % (tag, w["n"], g["n"])))
                    if not gd.same(g.get("r"), w["r"]):
                        bad.append((fi, "%s r: want %r got %r" % (tag, w["r"], g.get("r"))))
                    if not gd.same(g.get("statistic"), w["statistic"]):
                        bad.append((fi, "%s statistic: want %r got %r" % (tag, w["statistic"], g.get("statistic"))))
                    if g.get("df") != w["df"]:
                        bad.append((fi, "%s df: want %r got %r" % (tag, w["df"], g.get("df"))))
                    if w["p"] is None:
                        if g.get("p_value") is not None:
                            bad.append((fi, "%s p_value: want null got %r" % (tag, g.get("p_value"))))
                    elif not p_ok(g.get("p_value"), w["p"], w["n"]):
                        bad.append((fi, "%s p_value: want %s got %r" % (tag, mpmath.nstr(w["p"], 17), g.get("p_value"))))
                    if g.get("reason") != w["reason"]:
                        bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
                    if g.get("p_exact_zero", False) != w["one"]:
                        bad.append((fi, "%s p_exact_zero: want %r got %r" % (tag, w["one"], g.get("p_exact_zero"))))
                    if g.get("p_value") is not None and w["p"] is not None and w["p"] > mpmath.mpf("1e-300"):
                        err = float(abs(mpmath.mpf(g["p_value"]) - w["p"]) / w["p"])
                        pst["max"] = max(pst["max"], err)
                        pst["ratio"] = max(pst["ratio"], err / (1e-12 + (w["n"] - 1) * EPS))
                        pst["n"] += 1
                    if g.get("p_value") is not None and not (0.0 <= g["p_value"] <= 1.0):
                        bad.append((fi, "%s p_value %r is not in [0, 1]" % (tag, g["p_value"])))
                    if g.get("p_method") != "t":
                        bad.append((fi, "%s p_method %r" % (tag, g.get("p_method"))))
                    # the controls: a wrong degrees of freedom, a one-sided tail, a float one-pass r must each leave the bound on the cases where they differ
                    if w["p"] is not None and w["p"] > mpmath.mpf("1e-250") and not w["one"] and w["n"] >= 4:
                        nn_ = w["n"]
                        xq = Fraction(1) - Fraction(w["r"]) ** 2
                        wrong = mpmath.betainc(mpmath.mpf(nn_) / 2, mpmath.mpf(1) / 2, 0, mpmath.mpf(xq.numerator) / xq.denominator, regularized=True)
                        ctl["df-n"][1] += 1
                        ctl["df-n"][0] += 0 if p_ok(float(wrong), w["p"], nn_) else 1
                        ctl["one-sided"][1] += 1
                        ctl["one-sided"][0] += 0 if p_ok(float(w["p"] / 2), w["p"], nn_) else 1
                    if w["r"] is not None and w["n"] >= 3:
                        wr = welford_r(xs, ys)
                        ctl["welford"][1] += 1
                        ctl["welford"][0] += 0 if (wr is not None and wr == w["r"]) else 1
                    if a.scipy and w["r"] is not None and w["n"] >= 3 and not w["one"] and g.get("p_value") is not None:
                        import scipy.stats as ss
                        res = ss.pearsonr(xs, ys)
                        dr = abs(res[0] - g["r"])
                        dp = abs(res[1] - g["p_value"]) / max(abs(res[1]), 1e-300)
                        if shapes[x] in GOOD and shapes[y] in GOOD and abs(g["r"]) <= 0.999:
                            sp["n"] += 1
                            sp["r"] = max(sp["r"], dr)
                            sp["p"] = max(sp["p"], dp)
                            if dr > 1e-12 or dp > 1e-9:
                                bad.append((fi, "%s scipy: dr %g dp %g" % (tag, dr, dp)))
                        else:
                            sp["skipped"] += 1
                            sp["skipped_p"] = max(sp.get("skipped_p", 0.0), dp)
                # the order of the rows and of the columns
                if mode == "pairwise" and n > 1:
                    perm = list(range(n))
                    for i in range(n - 1, 0, -1):
                        k = r.next() % (i + 1)
                        perm[i], perm[k] = perm[k], perm[i]
                    p2 = run(a.bin, tmp, names, [rows[i] for i in perm], ["--missing", mode])
                    if p2.stdout != p.stdout:
                        bad.append((fi, "the rows shuffled give different bytes"))
                    # --cols names columns, the answer keeps the order of the file: the same bytes whatever the order given
                    sel = [names[j] for j in num]
                    p3 = run(a.bin, tmp, names, rows, ["--missing", mode, "--cols", ",".join(reversed(sel))])
                    if p3.stdout != p.stdout:
                        bad.append((fi, "--cols in another order gives different bytes"))
                    # the columns of the file in the opposite order: each pair is the mirrored one with the same values (r(x, y) = r(y, x) byte for byte)
                    names_r = list(reversed(names))
                    rows_r = [list(reversed(rw)) for rw in rows]
                    p4 = run(a.bin, tmp, names_r, rows_r, ["--missing", mode])
                    if p4.returncode != 0:
                        bad.append((fi, "mirrored file: exit %d" % p4.returncode))
                    else:
                        by_name = {frozenset((g["x"], g["y"])): g for g in got_pairs}
                        for g4 in json.loads(p4.stdout)["data"]["pairs"]:
                            g0 = by_name.get(frozenset((g4["x"], g4["y"])))
                            for k in ("n", "r", "statistic", "df", "p_value", "reason"):
                                if g0 is None or not gd.same(g0.get(k), g4.get(k)):
                                    bad.append((fi, "mirrored file: %s~%s key %s differs" % (g4["x"], g4["y"], k)))
                                    break
    for k, (outside, total) in ctl.items():
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_corr: %d files, %d pairs compared, %d differ; p-values against mpmath: largest relative error %.3g over %d, at most %.3g of the bound; controls outside the bound: %s%s" % (
        a.files, checked, len(bad), pst["max"], pst["n"], pst["ratio"], ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in ctl.items()),
        ("; scipy on %d well-conditioned pairs: max |dr| %.3g, max rel dp %.3g (%d others reported only, max rel dp %.3g)" % (sp["n"], sp["r"], sp["p"], sp["skipped"], sp.get("skipped_p", 0.0))) if a.scipy else ""))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
