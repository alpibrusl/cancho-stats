#!/usr/bin/env python3
"""S3.2 gates T-A, T-B, T-C, T-D and T-H (docs/test.md section 5) for `stats test t2` (Welch and pooled) and `stats test paired`:

    python3 scripts/gate_test2.py [--bin build/stats] [--kind t2|paired] [--files 200] [--seed 1] [--scipy]

Generated CSV files (the shapes of gate_describe.py, a group column of two levels whose sizes go down to 0, 1 and 2, missing cells in the values and in the group) against scripts/oracle_test.py, which works from
`Fraction` and `mpmath` and shares no code with the program. For each test: n of each group (or the pairs used), the statistic, the degrees of freedom (Welch: the rational of Satterthwaite rounded once, pooled: n1 + n2 - 2),
Cohen's d, the difference of means and the variance ratio (bytes), the p-value (relative 1e-11 + (df + 1) 2^-53 against mpmath at the EXACT rational statistic and degrees of freedom), the endpoints of the interval
(absolute 1e-12 max(|estimate|, half-width), at least two units of 2^-1074), the reasons. Also: the rows shuffled give the same bytes; the columns named in another order give the same bytes; the two groups
swapped (the labels renamed so that the byte order of the levels flips) negate the statistic, the difference and the effect size and leave df, p and the width of the interval (for `paired`, the pair named the other way round).
CONTROLS from the same data must leave the bounds on at least a tenth of the cases, so the gate is able to fail: for Welch the pooled test, for pooled the Welch one, the degrees of freedom plus one, the degrees
of freedom truncated (Welch), a one-sided p-value, the normal for Student's t, for the interval the normal quantile and the degrees of freedom plus one, and for `paired` the two columns treated as independent samples.
With --scipy: `ttest_ind` (equal_var False and True) and `ttest_rel` agree to 1e-12 (relative above |t| = 1, absolute below) on the statistic, 1e-12 on the degrees of freedom and 1e-9 relative on the p-value on the
well-conditioned columns (uniform, ints); the other families are reported only.

Exit 0: every number is the expected one. Exit 1: a difference. Exit 2: the harness failed.
"""
import argparse
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
import gate_describe as gd  # noqa: E402
import oracle_test as ot  # noqa: E402
from gate_test import GOOD, LEVELS, SHAPES, numeric, parse  # noqa: E402

LABELS = [("a", "b"), ("ctl", "trt"), ("0", "1"), ("10", "9"), ("B", "a"), ("x y", "x z")]

# hand cases of t2: (groups of values for column x as (a, b), method) -- groups of size 1 and 2, a constant group, both constant, a huge pair, subnormals, equal sizes and variances (Welch = pooled)
FIXED2 = [
    ([1.5], [2.5, 3.5, 4.0]),
    ([1.0, 2.0], [3.0, 5.0]),
    ([3.0, 3.0, 3.0], [1.0, 2.0, 4.0, 8.0]),
    ([3.0, 3.0, 3.0], [5.0, 5.0]),
    ([1e300, 2e300, -3e300, 5e299], [1e300, 3e300, 4e300]),
    ([5e-324, 1e-323, 1.5e-323, 2e-323], [0.0, 5e-324, 1e-323]),
    ([1.0, 2.0, 3.0, 4.0, 5.0], [2.0, 3.0, 4.0, 5.0, 6.0]),
    ([1e8, 100000001.0, 100000002.0, 100000004.0], [100000003.0, 100000007.0, 100000011.0]),
    ([], [1.0, 2.0, 3.0]),
]

# hand cases of paired: lists of (a, b) with None for a missing cell
FIXED_P = [
    [(1.5, 2.5)],
    [(1.5, 2.5), (3.0, 1.0)],
    [(3.0, 3.0), (4.0, 4.0), (5.0, 5.0)],
    [(1.0, 2.0), (2.0, 4.0), (3.0, 5.0), (4.0, 9.0), (5.0, 9.5)],
    [(1e300, -1e300), (2e300, 5e299), (-3e300, 1e300)],
    [(5e-324, 1e-323), (1e-323, 3e-323), (1.5e-323, 0.0), (2e-323, 5e-324)],
    [(1e8, 100000001.0), (100000001.0, 100000003.0), (100000002.0, 100000002.0), (100000004.0, 100000009.0)],
    [(1.0, None), (None, 2.0), (3.0, 4.0), (4.0, 7.0), (5.0, 5.5)],
]


def fmt(v):
    return "" if v is None else repr(v)


def run(binary, tmp, kind, names, rows, extra):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "test", kind, "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=300)


def shuffled(r, seq):
    perm = list(range(len(seq)))
    for i in range(len(perm) - 1, 0, -1):
        k = r.next() % (i + 1)
        perm[i], perm[k] = perm[k], perm[i]
    return [seq[i] for i in perm]


class Stat:
    def __init__(self):
        self.ctl = {}
        self.p = {"max": 0.0, "ratio": 0.0, "n": 0, "iv": 0.0, "iv_at": ""}
        self.sp = {"t": 0.0, "p": 0.0, "df": 0.0, "n": 0, "skipped": 0}
        self.checked = 0

    def control(self, name, outside):
        c = self.ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0


def check_common(bad, fi, tag, g, w, level, st, df_for_bound, header_of):
    """the numbers every kind answers: statistic, p-value, effect size, interval, reason"""
    if not gd.same(g.get("statistic"), w["statistic"]):
        bad.append((fi, "%s statistic: want %r got %r" % (tag, w["statistic"], g.get("statistic"))))
    es = g.get("effect_size")
    if w["d"] is None:
        if es is not None:
            bad.append((fi, "%s effect_size: want null got %r" % (tag, es)))
    elif es is None or es.get("name") != "cohen_d" or not gd.same(es.get("value"), w["d"]):
        bad.append((fi, "%s effect_size: want %r got %r" % (tag, w["d"], es)))
    if w["p"] is None:
        if g.get("p_value") is not None:
            bad.append((fi, "%s p_value: want null got %r" % (tag, g.get("p_value"))))
    elif not ot.p_ok(g.get("p_value"), w["p"], df_for_bound):
        bad.append((fi, "%s p_value: want %s got %r" % (tag, mpmath.nstr(w["p"], 17), g.get("p_value"))))
    elif w["p"] > mpmath.mpf("1e-300"):
        err = float(abs(mpmath.mpf(g["p_value"]) - w["p"]) / w["p"])
        st.p["max"] = max(st.p["max"], err)
        st.p["ratio"] = max(st.p["ratio"], err / (1e-11 + (df_for_bound + 1) * ot.EPS))
        st.p["n"] += 1
    if g.get("p_value") is not None and not (0.0 <= g["p_value"] <= 1.0):
        bad.append((fi, "%s p_value %r is not in [0, 1]" % (tag, g["p_value"])))
    if g.get("reason") != w["reason"]:
        bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
    if g.get("p_method") != "t":
        bad.append((fi, "%s p_method %r" % (tag, g.get("p_method"))))
    iv = g.get("interval")
    if w["low"] is None:
        if iv is not None:
            bad.append((fi, "%s interval: want null got %r" % (tag, iv)))
        return
    want_level = float(Fraction(ot.level_parts(level)[0], 10 ** ot.level_parts(level)[1]))
    if iv is None or iv.get("of") != header_of or iv.get("rounding") != "double" or not gd.same(iv.get("level"), want_level) or not gd.same(iv.get("estimate"), w["estimate"]):
        bad.append((fi, "%s interval header: want of %r estimate %r level %r got %r" % (tag, header_of, w["estimate"], want_level, iv)))
        return
    for k in ("low", "high"):
        if not ot.endpoint_ok(iv.get(k), w[k], w["estimate"], w["half"], w["q"]):
            bad.append((fi, "%s interval %s: want %s got %r" % (tag, k, mpmath.nstr(w[k], 17), iv.get(k))))
        else:
            scale = max(abs(mpmath.mpf(w["estimate"])), abs(w["half"]))
            ratio_iv = float(abs(mpmath.mpf(iv[k]) - w[k]) / scale) if scale > mpmath.mpf("1e-290") else 0.0
            if ratio_iv > st.p["iv"]:
                st.p["iv"] = ratio_iv
                st.p["iv_at"] = "%s %s est=%r half=%s got=%r want=%s" % (tag, k, w["estimate"], mpmath.nstr(w["half"], 6), iv[k], mpmath.nstr(w[k], 17))


def interval_controls(st, g, w, level, df_exact):
    """the interval of the answer against the endpoints of a wrong quantile: the normal one and the one for df + 1; only where the width is positive"""
    if w["low"] is None or not w["half"] > 0 or g.get("interval") is None:
        return
    zq = mpmath.sqrt(2) * mpmath.erfinv(1 - 2 * ot.frac_mpf(ot.alpha_half(level)))
    se = w["half"] / w["q"]   # the exact standard error the oracle's half-width is made of
    st.control("interval-normal-quantile", not ot.endpoint_ok(g["interval"]["high"], mpmath.mpf(w["estimate"]) + zq * se, w["estimate"], w["half"], w["q"]))
    qn = ot.t_quantile_upper(ot.alpha_half(level), ot.frac_mpf(df_exact) + 1)
    st.control("interval-df-plus-one", not ot.endpoint_ok(g["interval"]["high"], mpmath.mpf(w["estimate"]) + qn * se, w["estimate"], w["half"], w["q"]))


def p_controls(st, g, w, df_exact, pooled_p=None):
    if w["p"] is None or not w["p"] > mpmath.mpf("1e-250") or w["t2"] == 0:
        return
    dfm = ot.frac_mpf(df_exact)
    bound_df = float(df_exact)
    t2 = w["t2"]
    wrong = ot.p_two_sided_t(dfm + 1, Fraction(df_exact + 1) / (df_exact + 1 + t2))
    st.control("df-plus-one", not ot.p_ok(float(wrong), w["p"], bound_df))
    st.control("one-sided", not ot.p_ok(float(w["p"] / 2), w["p"], bound_df))
    tt = mpmath.sqrt(ot.frac_mpf(t2))
    st.control("normal", not ot.p_ok(float(mpmath.erfc(tt / mpmath.sqrt(2))), w["p"], bound_df))
    if df_exact.denominator != 1:
        dtr = Fraction(math.floor(df_exact))
        wrong = ot.p_two_sided_t(ot.frac_mpf(dtr), dtr / (dtr + t2))
        st.control("df-truncated", not ot.p_ok(float(wrong), w["p"], bound_df))
    if pooled_p is not None:
        st.control("other-variance-assumption", not ot.p_ok(float(pooled_p), w["p"], bound_df))


def scipy_check(bad, fi, tag, st, shape, g, res, welch_df=None):
    dt = abs(res.statistic - g["statistic"]) / max(abs(res.statistic), 1.0)
    dp = abs(res.pvalue - g["p_value"]) / max(abs(res.pvalue), 1e-300)
    ddf = abs(float(res.df) - float(g["df"])) / max(abs(float(res.df)), 1.0)
    if shape in GOOD:
        st.sp["n"] += 1
        st.sp["t"] = max(st.sp["t"], dt)
        st.sp["p"] = max(st.sp["p"], dp)
        st.sp["df"] = max(st.sp["df"], ddf)
        if dt > 1e-12 or dp > 1e-9 or ddf > 1e-12:
            bad.append((fi, "%s scipy: dt %g dp %g ddf %g" % (tag, dt, dp, ddf)))
    else:
        st.sp["skipped"] += 1


def groups_of(rows, gcol, labels):
    a = [i for i, rw in enumerate(rows) if rw[gcol] == labels[0]]
    b = [i for i, rw in enumerate(rows) if rw[gcol] == labels[1]]
    return a, b


def gen_t2(r, fi):
    """names, rows (strings), the group column index, the labels, the value columns' shapes and indices"""
    if fi < len(FIXED2):
        ga, gb = FIXED2[fi]
        labels = LABELS[fi % len(LABELS)]
        names = ["x", "g"]
        rows = [[fmt(v), labels[0]] for v in ga] + [[fmt(v), labels[1]] for v in gb]
        shapes = ["fixed", "group"]
        return names, rows, 1, sorted(labels), shapes, [0], r.pick(["welch", "pooled"]), r.pick(LEVELS)
    labels = sorted(LABELS[int(r.next() % len(LABELS))])
    nval = 1 + int(r.next() % 3)
    shapes = [SHAPES[int(r.next() % len(SHAPES))] for _ in range(nval)]
    n = [2, 3, 4, 5, 8, 20, 60, 200][int(r.next() % 8)]
    cols = [gd.column(r, s, n) for s in shapes]
    pa = [0.5, 0.5, 0.2, 0.8, 0.08][int(r.next() % 5)]
    grp = []
    for i in range(n):
        u = (r.next() % 10000) / 10000.0
        if u < 0.03:
            grp.append("")
        else:
            grp.append(labels[0] if (r.next() % 10000) / 10000.0 < pa else labels[1])
    if n >= 2:
        grp[0], grp[1] = labels[0], labels[1]
    gcol = int(r.next() % (nval + 1))
    names = ["c%d_%s" % (j, s) for j, s in enumerate(shapes)]
    allcols = cols[:]
    names.insert(gcol, "g")
    allcols.insert(gcol, grp)
    rows = [[allcols[j][i] for j in range(len(allcols))] for i in range(n)]
    vidx = [j for j in range(len(allcols)) if j != gcol]
    sh = shapes[:]
    sh.insert(gcol, "group")
    return names, rows, gcol, labels, sh, vidx, r.pick(["welch", "pooled"]), r.pick(LEVELS)


def test_t2(a, r, bad, st, tmp):
    files = 0
    for fi in range(a.files):
        names, rows, gcol, labels, shapes, vidx, method, level = gen_t2(r, fi)
        n = len(rows)
        num = [j for j in vidx if numeric([rw[j] for rw in rows])]
        if not num:
            continue
        chosen = [j for j in num if r.next() % 3 != 0] or [num[0]]
        pooled = method == "pooled"
        extra = ["--cols", ",".join(names[j] for j in chosen), "--group", "g", "--variance", method, "--level", level]
        ia, ib = groups_of(rows, gcol, labels)
        if fi < len(FIXED2) and not (ia and ib):
            # a group of no rows is a refusal (a group column with one level): covered by the refusals gate
            continue
        try:
            p = run(a.bin, tmp, "t2", names, rows, extra)
        except subprocess.TimeoutExpired:
            print("gate_test2: the program did not finish (hang)")
            return None
        if p.returncode != 0:
            keep = pathlib.Path(tempfile.gettempdir()) / ("gate_test2_failed_%d.csv" % fi)
            keep.write_text((tmp / "f.csv").read_text())
            print("gate_test2: exit %d on file %d: %s (the file is %s, args %r)" % (p.returncode, fi, p.stdout[:240], keep, extra))
            return None
        files += 1
        doc = json.loads(p.stdout)["data"]
        want_level = float(Fraction(ot.level_parts(level)[0], 10 ** ot.level_parts(level)[1]))
        if doc.get("test") != "t2" or doc.get("rows") != n or doc.get("group") != "g" or doc.get("levels") != {"a": labels[0], "b": labels[1]} or not gd.same(doc.get("level"), want_level) or doc.get("variance") != method or doc.get("adjustment") != "none":
            bad.append((fi, "the header: want t2 rows %d group g levels %r variance %s level %r; got %r" % (n, labels, method, want_level, {k: doc.get(k) for k in ("test", "rows", "group", "levels", "variance", "level", "adjustment")})))
        tests = doc["tests"]
        if [x["column"] for x in tests] != [names[j] for j in chosen]:
            bad.append((fi, "the tests are %r" % [x["column"] for x in tests]))
            continue
        for g, j in zip(tests, chosen):
            va = [v for v in (parse(rows[i][j]) for i in ia) if v is not None]
            vb = [v for v in (parse(rows[i][j]) for i in ib) if v is not None]
            w = ot.t2(va, vb, pooled, level)
            other = ot.t2(va, vb, not pooled, level)
            tag = "t2 %s %s level %s" % (names[j], method, level)
            st.checked += 1
            if g["n"] != {"a": w["n1"], "b": w["n2"]}:
                bad.append((fi, "%s n: want %r got %r" % (tag, {"a": w["n1"], "b": w["n2"]}, g["n"])))
            if g.get("method") != ("pooled" if pooled else "welch"):
                bad.append((fi, "%s method %r" % (tag, g.get("method"))))
            if w["df"] is None:
                if g.get("df") is not None:
                    bad.append((fi, "%s df: want null got %r" % (tag, g.get("df"))))
            elif pooled:
                if g.get("df") != int(w["df"]) or not isinstance(g.get("df"), int):
                    bad.append((fi, "%s df: want %d got %r" % (tag, w["df"], g.get("df"))))
            elif not gd.same(g.get("df"), ot.rounded(w["df"])):
                bad.append((fi, "%s df: want %r got %r" % (tag, ot.rounded(w["df"]), g.get("df"))))
            vr = (g.get("assumptions") or {}).get("variance_ratio", "missing")
            if (w["variance_ratio"] is None and vr is not None) or (w["variance_ratio"] is not None and not gd.same(vr, w["variance_ratio"])):
                bad.append((fi, "%s variance_ratio: want %r got %r" % (tag, w["variance_ratio"], vr)))
            bound_df = float(w["df"]) if w["df"] is not None else 1.0
            check_common(bad, fi, tag, g, w, level, st, bound_df, "difference of means")
            if w["p"] is not None:
                p_controls(st, g, w, w["df"], other["p"] if other["p"] is not None else None)
                interval_controls(st, g, w, level, w["df"])
            elif w["low"] is not None and w["half"] > 0 and w["df"] is not None:
                interval_controls(st, g, w, level, w["df"])
            if a.scipy and w["statistic"] is not None and w["p"] is not None and g.get("p_value") is not None:
                import scipy.stats as ss
                res = ss.ttest_ind(va, vb, equal_var=pooled)
                scipy_check(bad, fi, tag, st, shapes[j], g, res)
        # invariance
        if n > 1:
            p2 = run(a.bin, tmp, "t2", names, shuffled(r, rows), extra)
            if p2.stdout != p.stdout:
                bad.append((fi, "the rows shuffled give different bytes"))
            if len(chosen) > 1:
                p3 = run(a.bin, tmp, "t2", names, rows, ["--cols", ",".join(reversed([names[j] for j in chosen]))] + extra[2:])
                if p3.stdout != p.stdout:
                    bad.append((fi, "--cols in another order gives different bytes"))
            # the groups swapped: the labels renamed so that the byte order flips
            swap = {labels[0]: "zz" + labels[1], labels[1]: "yy" + labels[0]}
            srows = [[(swap[c] if k == gcol and c in swap else c) for k, c in enumerate(rw)] for rw in rows]
            p4 = run(a.bin, tmp, "t2", names, srows, extra)
            if p4.returncode != 0:
                bad.append((fi, "swapped groups: exit %d" % p4.returncode))
            else:
                d4 = json.loads(p4.stdout)["data"]["tests"]
                for g0, g4 in zip(tests, d4):
                    ok = g0["df"] == g4["df"] and g0["p_value"] == g4["p_value"] and g0.get("reason") == g4.get("reason") and g0["n"] == {"a": g4["n"]["b"], "b": g4["n"]["a"]}
                    for k in ("statistic",):
                        if g0[k] is not None:
                            ok = ok and gd.same(-g0[k] if g0[k] != 0 else 0.0, g4[k])
                    if g0.get("effect_size") is not None:
                        v0, v4 = g0["effect_size"]["value"], g4["effect_size"]["value"]
                        ok = ok and gd.same(-v0 if v0 != 0 else 0.0, v4)
                    iv0, iv4 = g0.get("interval"), g4.get("interval")
                    if iv0 is not None:
                        ok = ok and gd.same(-iv0["estimate"] if iv0["estimate"] != 0 else 0.0, iv4["estimate"]) and gd.same(-iv0["high"] if iv0["high"] != 0 else 0.0, iv4["low"]) and gd.same(-iv0["low"] if iv0["low"] != 0 else 0.0, iv4["high"])
                    if not ok:
                        bad.append((fi, "the groups swapped do not negate the answer for %s: %r against %r" % (g0["column"], g0, g4)))
    return files


def gen_paired(r, fi):
    if fi < len(FIXED_P):
        pr = FIXED_P[fi]
        names = ["a", "b"]
        rows = [[fmt(x), fmt(y)] for x, y in pr]
        return names, rows, [("a", "b")], ["fixed", "fixed"], r.pick(LEVELS)
    ncols = 2 + int(r.next() % 3)
    shapes = [SHAPES[int(r.next() % len(SHAPES))] for _ in range(ncols)]
    n = [2, 3, 4, 5, 8, 20, 60, 200][int(r.next() % 8)]
    cols = [gd.column(r, s, n) for s in shapes]
    names = ["c%d_%s" % (j, s) for j, s in enumerate(shapes)]
    rows = [[cols[j][i] for j in range(ncols)] for i in range(n)]
    num = [j for j in range(ncols) if numeric(cols[j])]
    pairs = []
    for _ in range(1 + int(r.next() % 3)):
        if len(num) >= 2:
            x = num[int(r.next() % len(num))]
            y = num[int(r.next() % len(num))]
            if x == y and r.next() % 4 != 0:
                y = num[(num.index(x) + 1) % len(num)]
            pairs.append((names[x], names[y]))
    return names, rows, pairs, shapes, r.pick(LEVELS)


def test_paired(a, r, bad, st, tmp):
    files = 0
    for fi in range(a.files):
        names, rows, pairs, shapes, level = gen_paired(r, fi)
        if not pairs:
            continue
        n = len(rows)
        seen = []
        for pr in pairs:
            if pr not in seen:
                seen.append(pr)
        pairs = seen
        extra = ["--pairs", ",".join("%s:%s" % pr for pr in pairs), "--level", level]
        p = run(a.bin, tmp, "paired", names, rows, extra)
        if p.returncode != 0:
            keep = pathlib.Path(tempfile.gettempdir()) / ("gate_test2_failed_p%d.csv" % fi)
            keep.write_text((tmp / "f.csv").read_text())
            print("gate_test2: exit %d on paired file %d: %s (the file is %s, args %r)" % (p.returncode, fi, p.stdout[:240], keep, extra))
            return None
        files += 1
        doc = json.loads(p.stdout)["data"]
        want_level = float(Fraction(ot.level_parts(level)[0], 10 ** ot.level_parts(level)[1]))
        if doc.get("test") != "paired" or doc.get("rows") != n or not gd.same(doc.get("level"), want_level) or doc.get("adjustment") != "none":
            bad.append((fi, "the header: want paired rows %d level %r; got %r" % (n, want_level, {k: doc.get(k) for k in ("test", "rows", "level", "adjustment")})))
        tests = doc["tests"]
        if [x["pair"] for x in tests] != [list(pr) for pr in pairs]:
            bad.append((fi, "the tests are %r" % [x["pair"] for x in tests]))
            continue
        for g, (na, nb) in zip(tests, pairs):
            ja, jb = names.index(na), names.index(nb)
            used = [(parse(rw[ja]), parse(rw[jb])) for rw in rows]
            used = [(x, y) for x, y in used if x is not None and y is not None]
            w = ot.paired(used, level)
            tag = "paired %s:%s level %s" % (na, nb, level)
            st.checked += 1
            if g["n"] != {"n": w["n"], "missing": n - w["n"]}:
                bad.append((fi, "%s n: want %r got %r" % (tag, {"n": w["n"], "missing": n - w["n"]}, g["n"])))
            if g.get("df") != w["df"]:
                bad.append((fi, "%s df: want %r got %r" % (tag, w["df"], g.get("df"))))
            if g.get("method") != "paired t":
                bad.append((fi, "%s method %r" % (tag, g.get("method"))))
            bound_df = w["df"] if w["df"] is not None else 1
            check_common(bad, fi, tag, g, w, level, st, bound_df, "difference of means")
            if w["p"] is not None and w["n"] >= 3:
                p_controls(st, g, w, Fraction(w["df"]))
                # the two columns treated as independent samples
                xs = [x for x, _ in used]
                ys = [y for _, y in used]
                ind = ot.t2(xs, ys, False, level)
                if ind["p"] is not None:
                    st.control("independent-samples", not ot.p_ok(float(ind["p"]), w["p"], bound_df))
            if w["low"] is not None and w["half"] > 0 and w["n"] >= 3:
                interval_controls(st, g, w, level, Fraction(w["df"]))
            if a.scipy and w["statistic"] is not None and w["p"] is not None and g.get("p_value") is not None:
                import scipy.stats as ss
                res = ss.ttest_rel([x for x, _ in used], [y for _, y in used])
                scipy_check(bad, fi, tag, st, shapes[ja] if shapes[ja] == shapes[jb] else "mixed", g, res)
        if n > 1:
            p2 = run(a.bin, tmp, "paired", names, shuffled(r, rows), extra)
            if p2.stdout != p.stdout:
                bad.append((fi, "the rows shuffled give different bytes"))
            p3 = run(a.bin, tmp, "paired", names, rows, ["--pairs", ",".join("%s:%s" % pr for pr in reversed(pairs)), "--level", level])
            if len(pairs) > 1 and [x["pair"] for x in json.loads(p3.stdout)["data"]["tests"]] != [list(pr) for pr in reversed(pairs)]:
                bad.append((fi, "--pairs in another order does not keep the order named"))
            # the pairs named the other way round negate the answer
            flipped = ["%s:%s" % (pb, pa) for pa, pb in pairs]
            p4 = run(a.bin, tmp, "paired", names, rows, ["--pairs", ",".join(flipped), "--level", level])
            d4 = json.loads(p4.stdout)["data"]["tests"] if p4.returncode == 0 else None
            if d4 is None:
                bad.append((fi, "flipped pairs: exit %d" % p4.returncode))
            else:
                for g0, g4 in zip(tests, d4):
                    ok = g0["df"] == g4["df"] and g0["p_value"] == g4["p_value"] and g0.get("reason") == g4.get("reason") and g0["n"] == g4["n"]
                    if g0["statistic"] is not None:
                        ok = ok and gd.same(-g0["statistic"] if g0["statistic"] != 0 else 0.0, g4["statistic"])
                    if g0.get("effect_size") is not None:
                        v0, v4 = g0["effect_size"]["value"], g4["effect_size"]["value"]
                        ok = ok and gd.same(-v0 if v0 != 0 else 0.0, v4)
                    iv0, iv4 = g0.get("interval"), g4.get("interval")
                    if iv0 is not None:
                        ok = ok and gd.same(-iv0["estimate"] if iv0["estimate"] != 0 else 0.0, iv4["estimate"]) and gd.same(-iv0["high"] if iv0["high"] != 0 else 0.0, iv4["low"])
                    if not ok:
                        bad.append((fi, "the flipped pair does not negate the answer: %r against %r" % (g0, g4)))
    return files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--kind", default="t2", choices=["t2", "paired"])
    ap.add_argument("--files", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scipy", action="store_true")
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    st = Stat()
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        files = (test_t2 if a.kind == "t2" else test_paired)(a, r, bad, st, tmp)
    if files is None:
        return 2
    for k, (outside, total) in sorted(st.ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_test2 %s: %d files, %d tests compared, %d differ; p-values against mpmath: largest relative error %.3g over %d, at most %.3g of the bound; interval endpoints: largest error %.3g of max(|estimate|, half-width) (bound 1e-12); controls outside the bound: %s%s" % (
        a.kind, a.files, st.checked, len(bad), st.p["max"], st.p["n"], st.p["ratio"], st.p["iv"], ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(st.ctl.items())),
        ("; scipy on %d well-conditioned columns: max relative dt %.3g, ddf %.3g, dp %.3g (%d others not gated)" % (st.sp["n"], st.sp["t"], st.sp["df"], st.sp["p"], st.sp["skipped"])) if a.scipy else ""))
    if st.p["iv_at"]:
        print("  largest interval error at: %s" % st.p["iv_at"])
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
