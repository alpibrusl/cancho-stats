#!/usr/bin/env python3
"""S2.2 and S2.3 gates C-B, C-C and C-D (docs/corr.md section 5), for `stats corr` (Pearson, with --method spearman the average-rank coefficient, with --method kendall tau-b): generated CSV files against

    python3 scripts/gate_corr.py [--bin build/stats] [--files 200] [--seed 1] [--scipy] [--method pearson|spearman|kendall]

With --method spearman the oracle ranks the rows used for each pair by counting (twice the average rank = 2 (values below) + (values equal, itself included) + 1, an independent method from the program's
binary search), takes the exact Pearson coefficient of those ranks, and expects the key `rho`, the p-method `t-approximation`, the tie groups per column, and the reason `all-tied` for a column of equal ranks. With --method kendall the oracle counts concordant, discordant and tied pairs over ALL pairs (O(n^2)), takes the tie-corrected variance from the group sizes in Kendall's textbook form, and expects the key `tau`,
`s` (= C - D), `ties` with the rule `tau-b`, the p-method `normal-approximation`, no degrees of freedom, and the reasons `all-tied`, `too-few-rows` (n < 3: no statistic) and `variance-zero`; its p-value bound is relative
1e-12 + (z^2 + 1) 2^-53 (the rounding of z amplified by z^2). Controls: the untied variance on tied data, tau-a, a one-sided tail. For spearman: ties broken
by position, and ranks taken over the whole column before the rows with a missing partner are dropped; each must give another coefficient on the cases where it applies.

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
GOOD = ("uniform", "ints", "lin", "tied", "tied3")
SHAPES = ["uniform", "ints", "tied", "tied3", "tied3", "lin", "offset", "offset-neg", "wide", "tiny", "huge", "const", "gappy", "empty", "text", "mixed", "lin", "lin", "ints", "uniform"]


def _t(*cols):
    return [[repr(float(v)) for v in c] for c in cols]


# hand cases that every run includes: |t| beyond the largest double (n = 3, the tiny third row decides the last digits), r exactly -1 and +1, n = 2, subnormals, a huge exact pair
FIXED = [
    _t([1.0, 2.0, 2.0 ** -1000], [2.0, 4.0, 2 * 2.0 ** -1000 * (1 + 2.0 ** -52)]),
    # |t| about 1e257 with df = 1: the tail, 6.1e-258, is representable (the Student t survival function answered 0 here until the correction of docs/numerics.md section 3.2)
    _t([1.0, 2.0, 2.0 ** -800], [2.0, 4.0, 2 * 2.0 ** -800 * (1 + 2.0 ** -52)]),
    _t([1, 2, 3, 4], [-2, -4, -6, -8]),
    _t([1, 2], [3, 5]),
    _t([1e300, 2e300, -3e300], [1e300, 2e300, -3e300]),
    _t([5e-324, 1e-323, 1.5e-323, 2e-323], [1e-323, 3e-323, 2e-323, 4e-323]),
    _t([1.0, 2.0, 3.0], [1.0, 2.0, 3.0000000000000004]),
    _t([1e8, 1e8 + 1, 1e8 + 2, 1e8 + 4], [3.0, 1.0, 2.0, 5.0], [7.0, 7.0, 7.0, 7.0]),
    # groups of ties of sizes 5, 3, 2, 1, 1 in x and 4, 4, 2, 2 in y: every term of the variance of tau-b (v0, vt, vu, v1, v2) is nonzero and none is a special case
    _t([1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5], [1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4], [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8]),
    _t([1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5], [2, 2, 2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 1, 2, 3], [5, 4, 3, 2, 1, 5, 4, 3, 2, 1, 5, 4, 3, 2, 1]),
    _t([1, 2, 3, 4, 5], [1, 2, 3, 2, 1], [1, 2, 3, 2, 1.0000000000000002]),  # r = 0 exactly, and r about 1e-17: the p-value is 1, never above it
]


def column(r, shape, n, prev):
    """the cells of one column; `lin` follows the previous column (a correlated pair, with noise); `tied` takes four values (many ties)"""
    if shape == "tied3":
        # three values, so groups of ties of many sizes, in both columns of a pair: the terms of Kendall's variance that need groups of size 3 and more
        return [r.pick(["1", "1", "1", "2", "2", "3"]) for _ in range(n)]
    if shape == "tied":
        return [r.pick(["1", "2", "2.5", "4"]) if r.next() % 9 else "" for _ in range(n)]
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


def ranks2(vs):
    """twice the average rank of each value (1-based), by counting: 2 * (values below) + (values equal, itself included) + 1; -0.0 equals 0.0"""
    return [2 * sum(1 for u in vs if u < v) + sum(1 for u in vs if u == v) + 1 for v in vs]


def tie_groups(vs):
    """the number of values that occur more than once"""
    seen = {}
    for v in vs:
        seen[v] = seen.get(v, 0) + 1
    return sum(1 for c in seen.values() if c > 1)


def ordinal2(vs):
    """control: ties broken by position instead of averaged, doubled"""
    order = sorted(range(len(vs)), key=lambda i: (vs[i], i))
    out = [0] * len(vs)
    for rank, i in enumerate(order):
        out[i] = 2 * (rank + 1)
    return out


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


def expected_kendall(xs, ys):
    """the oracle of Kendall's tau-b for one pair, by the O(n^2) count over all pairs (no merge sort, no sorting): concordant, discordant, pairs tied in x, pairs tied in y; the tie-corrected variance
    of S = C - D from the sizes of the groups of ties, in its textbook form (Kendall 1970, not the common denominator the program uses); tau = S / sqrt((n0 - n1)(n0 - n2)) and z = S / sqrt(V) each rounded
    once, the two-sided normal p-value from the exact z at 60 digits."""
    n = len(xs)
    c = d = tx = ty = 0
    for i in range(n):
        for j in range(i + 1, n):
            sx = (xs[i] > xs[j]) - (xs[i] < xs[j])
            sy = (ys[i] > ys[j]) - (ys[i] < ys[j])
            if sx == 0:
                tx += 1
            if sy == 0:
                ty += 1
            if sx * sy > 0:
                c += 1
            elif sx * sy < 0:
                d += 1
    S = c - d
    n0 = n * (n - 1) // 2
    out = {"n": n, "r": None, "statistic": None, "df": None, "p": None, "one": False, "reason": None, "s": S, "zz": 0.0, "n0": n0, "tx": tx, "ty": ty}
    dx, dy = n0 - tx, n0 - ty
    if dx == 0 or dy == 0:
        out["reason"] = "all-tied"
        return out
    tau = gd.sqrt_round(Fraction(S * S, dx * dy))
    out["r"] = -tau if S < 0 else tau
    if n < 3:
        out["reason"] = "too-few-rows"
        return out
    def groups(vs):
        cnt = {}
        for v in vs:
            cnt[v] = cnt.get(v, 0) + 1
        return [t for t in cnt.values() if t > 1]
    gx, gy = groups(xs), groups(ys)
    v0 = n * (n - 1) * (2 * n + 5)
    vt = sum(t * (t - 1) * (2 * t + 5) for t in gx)
    vu = sum(t * (t - 1) * (2 * t + 5) for t in gy)
    v1 = sum(t * (t - 1) for t in gx) * sum(t * (t - 1) for t in gy)
    v2 = sum(t * (t - 1) * (t - 2) for t in gx) * sum(t * (t - 1) * (t - 2) for t in gy)
    V = Fraction(v0 - vt - vu, 18) + Fraction(v1, 2 * n * (n - 1)) + Fraction(v2, 9 * n * (n - 1) * (n - 2))
    out["V"] = V
    if V <= 0:
        out["reason"] = "variance-zero"
        return out
    z = gd.sqrt_round(Fraction(S * S) / V)
    out["statistic"] = -z if S < 0 else z
    out["zz"] = z * z
    zexact = mpmath.mpf(abs(S)) / mpmath.sqrt(mpmath.mpf(V.numerator) / V.denominator)
    out["p"] = mpmath.erfc(zexact / mpmath.sqrt(2))
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


def p_ok(got, want, n, z2=None):
    """the bound of C-C; for the normal approximation of Kendall the rounding of z is amplified by z^2 instead of n - 1"""
    if want == 0:
        return got == 0
    if want < mpmath.mpf("1e-300"):
        return got is not None and got < 1e-290
    slack = (n - 1) * EPS if z2 is None else (z2 + 1.0) * EPS
    return got is not None and abs(mpmath.mpf(got) - want) <= (mpmath.mpf("1e-12") + slack) * want


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
    ap.add_argument("--method", choices=["pearson", "spearman", "kendall"], default="pearson")
    a = ap.parse_args()
    r = SM64(a.seed)
    spear = a.method == "spearman"
    kend = a.method == "kendall"
    key = "rho" if spear else ("tau" if kend else "r")
    pmeth = "t-approximation" if spear else ("normal-approximation" if kend else "t")
    bad = []
    checked = 0
    ctl = {"df-n": [0, 0], "one-sided": [0, 0], ("ordinal-rank" if spear else "welford"): [0, 0]}  # [outside the bound, cases]
    if spear:
        ctl["rank-before-deletion"] = [0, 0]
    if kend:
        ctl = {"one-sided": [0, 0], "untied-variance": [0, 0], "tau-a": [0, 0]}
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
                    p = run(a.bin, tmp, names, rows, ["--missing", mode, "--method", a.method])
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
                    if spear:
                        w = expected([float(v) for v in ranks2(xs)], [float(v) for v in ranks2(ys)])
                        if w["reason"] == "constant-column":
                            w["reason"] = "all-tied"
                        want_ties = {"x": tie_groups(xs), "y": tie_groups(ys), "rule": "average"}
                    elif kend:
                        w = expected_kendall(xs, ys)
                        want_ties = {"x": tie_groups(xs), "y": tie_groups(ys), "rule": "tau-b"}
                    else:
                        w = expected(xs, ys)
                    tag = "%s %s~%s" % (mode, names[x], names[y])
                    checked += 1
                    if g["n"] != w["n"]:
                        bad.append((fi, "%s n: want %r got %r" % (tag, w["n"], g["n"])))
                    if not gd.same(g.get(key), w["r"]):
                        bad.append((fi, "%s %s: want %r got %r" % (tag, key, w["r"], g.get(key))))
                    if kend and g.get("s") != w["s"]:
                        bad.append((fi, "%s s: want %r got %r" % (tag, w["s"], g.get("s"))))
                    if (spear or kend) and g.get("ties") != want_ties:
                        bad.append((fi, "%s ties: want %r got %r" % (tag, want_ties, g.get("ties"))))
                    if not gd.same(g.get("statistic"), w["statistic"]):
                        bad.append((fi, "%s statistic: want %r got %r" % (tag, w["statistic"], g.get("statistic"))))
                    if g.get("df") != w["df"]:
                        bad.append((fi, "%s df: want %r got %r" % (tag, w["df"], g.get("df"))))
                    if w["p"] is None:
                        if g.get("p_value") is not None:
                            bad.append((fi, "%s p_value: want null got %r" % (tag, g.get("p_value"))))
                    elif not p_ok(g.get("p_value"), w["p"], w["n"], w["zz"] if kend else None):
                        bad.append((fi, "%s p_value: want %s got %r" % (tag, mpmath.nstr(w["p"], 17), g.get("p_value"))))
                    if g.get("reason") != w["reason"]:
                        bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
                    if g.get("p_exact_zero", False) != w["one"]:
                        bad.append((fi, "%s p_exact_zero: want %r got %r" % (tag, w["one"], g.get("p_exact_zero"))))
                    if g.get("p_value") is not None and w["p"] is not None and w["p"] > mpmath.mpf("1e-300"):
                        err = float(abs(mpmath.mpf(g["p_value"]) - w["p"]) / w["p"])
                        pst["max"] = max(pst["max"], err)
                        pst["ratio"] = max(pst["ratio"], err / (1e-12 + ((w["zz"] + 1.0) if kend else (w["n"] - 1)) * EPS))
                        pst["n"] += 1
                    if g.get("p_value") is not None and not (0.0 <= g["p_value"] <= 1.0):
                        bad.append((fi, "%s p_value %r is not in [0, 1]" % (tag, g["p_value"])))
                    if g.get("p_method") != pmeth:
                        bad.append((fi, "%s p_method %r" % (tag, g.get("p_method"))))
                    # the controls: a wrong degrees of freedom, a one-sided tail, a float one-pass r must each leave the bound on the cases where they differ
                    if kend and w["p"] is not None and w["p"] > mpmath.mpf("1e-250"):
                        ctl["one-sided"][1] += 1
                        ctl["one-sided"][0] += 0 if p_ok(float(w["p"] / 2), w["p"], w["n"], w["zz"]) else 1
                        if w["tx"] or w["ty"]:
                            # control: the variance of the untied case n (n - 1) (2 n + 5) / 18 where there are ties
                            nn_ = w["n"]
                            v_untied = Fraction(nn_ * (nn_ - 1) * (2 * nn_ + 5), 18)
                            zu = mpmath.mpf(abs(w["s"])) / mpmath.sqrt(mpmath.mpf(v_untied.numerator) / v_untied.denominator)
                            ctl["untied-variance"][1] += 1
                            ctl["untied-variance"][0] += 0 if p_ok(float(mpmath.erfc(zu / mpmath.sqrt(2))), w["p"], nn_, w["zz"]) else 1
                            # control: tau-a, S / n0, in place of tau-b
                            ta = gd.sqrt_round(Fraction(w["s"] * w["s"], w["n0"] * w["n0"]))
                            ctl["tau-a"][1] += 1
                            ctl["tau-a"][0] += 0 if gd.same(-ta if w["s"] < 0 else ta, w["r"]) else 1
                    if not kend and w["p"] is not None and w["p"] > mpmath.mpf("1e-250") and not w["one"] and w["n"] >= 4:
                        nn_ = w["n"]
                        xq = Fraction(1) - Fraction(w["r"]) ** 2
                        wrong = mpmath.betainc(mpmath.mpf(nn_) / 2, mpmath.mpf(1) / 2, 0, mpmath.mpf(xq.numerator) / xq.denominator, regularized=True)
                        ctl["df-n"][1] += 1
                        ctl["df-n"][0] += 0 if p_ok(float(wrong), w["p"], nn_) else 1
                        ctl["one-sided"][1] += 1
                        ctl["one-sided"][0] += 0 if p_ok(float(w["p"] / 2), w["p"], nn_) else 1
                    if w["r"] is not None and w["n"] >= 3 and not spear and not kend:
                        wr = welford_r(xs, ys)
                        ctl["welford"][1] += 1
                        ctl["welford"][0] += 0 if (wr is not None and wr == w["r"]) else 1
                    if spear and w["r"] is not None and w["n"] >= 3 and (want_ties["x"] or want_ties["y"]):
                        # control: ties broken by position; its coefficient must differ from the average-rank one
                        wo = expected([float(v) for v in ordinal2(xs)], [float(v) for v in ordinal2(ys)])
                        ctl["ordinal-rank"][1] += 1
                        ctl["ordinal-rank"][0] += 0 if gd.same(wo["r"], w["r"]) else 1
                    if spear and mode == "pairwise" and w["r"] is not None and w["n"] >= 3 and len(used) < n:
                        # control: the ranks taken over the whole column BEFORE the rows with a missing partner are dropped
                        fx_all = [vals[x][i] for i in range(n) if vals[x][i] is not None]
                        fy_all = [vals[y][i] for i in range(n) if vals[y][i] is not None]
                        rxa, rya = ranks2(fx_all), ranks2(fy_all)
                        ix = {i: k for k, i in enumerate([i for i in range(n) if vals[x][i] is not None])}
                        iy = {i: k for k, i in enumerate([i for i in range(n) if vals[y][i] is not None])}
                        wd = expected([float(rxa[ix[i]]) for i in used], [float(rya[iy[i]]) for i in used])
                        ctl["rank-before-deletion"][1] += 1
                        ctl["rank-before-deletion"][0] += 0 if gd.same(wd["r"], w["r"]) else 1
                    if a.scipy and w["r"] is not None and w["n"] >= 3 and not w["one"] and g.get("p_value") is not None:
                        import scipy.stats as ss
                        res = ss.spearmanr(xs, ys) if spear else (ss.kendalltau(xs, ys, variant="b", method="asymptotic") if kend else ss.pearsonr(xs, ys))
                        dr = abs(res[0] - g[key])
                        dp = abs(res[1] - g["p_value"]) / max(abs(res[1]), 1e-300)
                        if shapes[x] in GOOD and shapes[y] in GOOD and abs(g[key]) <= 0.999:
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
                    p2 = run(a.bin, tmp, names, [rows[i] for i in perm], ["--missing", mode, "--method", a.method])
                    if p2.stdout != p.stdout:
                        bad.append((fi, "the rows shuffled give different bytes"))
                    # --cols names columns, the answer keeps the order of the file: the same bytes whatever the order given
                    sel = [names[j] for j in num]
                    p3 = run(a.bin, tmp, names, rows, ["--missing", mode, "--method", a.method, "--cols", ",".join(reversed(sel))])
                    if p3.stdout != p.stdout:
                        bad.append((fi, "--cols in another order gives different bytes"))
                    # the columns of the file in the opposite order: each pair is the mirrored one with the same values (r(x, y) = r(y, x) byte for byte)
                    names_r = list(reversed(names))
                    rows_r = [list(reversed(rw)) for rw in rows]
                    p4 = run(a.bin, tmp, names_r, rows_r, ["--missing", mode, "--method", a.method])
                    if p4.returncode != 0:
                        bad.append((fi, "mirrored file: exit %d" % p4.returncode))
                    else:
                        by_name = {frozenset((g["x"], g["y"])): g for g in got_pairs}
                        for g4 in json.loads(p4.stdout)["data"]["pairs"]:
                            g0 = by_name.get(frozenset((g4["x"], g4["y"])))
                            for k in ("n", key, "statistic", "df", "p_value", "reason"):
                                if g0 is None or not gd.same(g0.get(k), g4.get(k)):
                                    bad.append((fi, "mirrored file: %s~%s key %s differs" % (g4["x"], g4["y"], k)))
                                    break
    if kend:
        # the variance of S is positive whenever tau-b is defined and n >= 3: exhaustively for every pair of columns of n <= 5 rows over the values 1, 2, 3 (the branch `variance-zero` of the program is then
        # a refusal that no input of this range reaches); checked on the oracle, and a sample of those inputs through the program
        import itertools
        count = bad_v = 0
        sample = []
        for nn_ in range(3, 6):
            for xs_ in itertools.product((1.0, 2.0, 3.0), repeat=nn_):
                for ys_ in itertools.product((1.0, 2.0, 3.0), repeat=nn_):
                    w_ = expected_kendall(list(xs_), list(ys_))
                    if w_["reason"] is None or w_["reason"] == "variance-zero":
                        count += 1
                        if w_["reason"] == "variance-zero":
                            bad_v += 1
                            bad.append((-1, "variance-zero for x=%r y=%r" % (xs_, ys_)))
                    if r.next() % 997 == 0:
                        sample.append((xs_, ys_))
        with tempfile.TemporaryDirectory() as t2:
            tmp2 = pathlib.Path(t2)
            for xs_, ys_ in sample[:200]:
                w_ = expected_kendall(list(xs_), list(ys_))
                p_ = run(a.bin, tmp2, ["x", "y"], [[repr(u), repr(v)] for u, v in zip(xs_, ys_)], ["--method", "kendall"])
                g_ = json.loads(p_.stdout)["data"]["pairs"][0]
                checked += 1
                if g_.get("reason") != w_["reason"] or not gd.same(g_.get("tau"), w_["r"]) or not gd.same(g_.get("statistic"), w_["statistic"]) or g_.get("s") != w_["s"]:
                    bad.append((-1, "exhaustive sample x=%r y=%r: want %r got %r" % (xs_, ys_, (w_["reason"], w_["r"], w_["statistic"], w_["s"]), g_)))
        print("gate_corr: variance of S positive for all %d pairs of columns with tau-b defined (n 3 to 5, values 1 to 3): %d violations; %d of them run through the program" % (count, bad_v, min(len(sample), 200)))
    for k, (outside, total) in ctl.items():
        if k in ("ordinal-rank", "rank-before-deletion", "untied-variance", "tau-a") and total < 5:
            bad.append((-1, "control %s had only %d applicable cases" % (k, total)))
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
