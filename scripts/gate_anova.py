#!/usr/bin/env python3
"""S4.1 gates R-A, R-B and R-G for `stats test anova` (docs/rank.md sections 3 and 4):

    python3 scripts/gate_anova.py [--bin build/stats] [--files 200] [--seed 1] [--scipy]

Generated CSV files (the shapes of gate_describe.py, a group column of 2 to 25 levels whose sizes go down to 0, 1 and 2, missing cells in the values and in the group, groups with their own mean and spread) against
scripts/oracle_test.py, which works from `Fraction` and `mpmath` and shares no code with the program. For each test and for both variance assumptions: n, mean and sd of every group (bytes), the statistic (the classic F is one
rounding of an exact rational: bytes; Welch's F* is double-double: bytes too, and the cases where it is one unit off are counted and reported), the degrees of freedom, eta^2 and omega^2 (bytes), the p-value (relative
1e-11 + (d1 + d2 + 1) 2^-53 against mpmath at the EXACT rational statistic and degrees of freedom), the reasons. Also: the rows shuffled give the same bytes; the columns named in another order give the same bytes; the levels
renamed so that their byte order is reversed give the same statistic, p-value and effect sizes (Welch's within one unit) and the per-group lists reversed.
CONTROLS from the same data must leave the bounds on at least a tenth of the cases, so the gate is able to fail: the degrees of freedom swapped, d2 + 1, the other variance assumption, and the chi-square approximation of F d1.
With --scipy: `f_oneway` agrees to 1e-10 + 1e-14 F relative on the statistic and 1e-8 on the p-value (scipy sums floats: its own error is 3e-12 in general and 2e-10 at F = 7e4, where the within-group sum of squares is a small difference) on the well-conditioned columns (uniform, ints, groupmix) for the classic F.

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
from gate_test import GOOD, parse  # noqa: E402

SHAPES = ["groupmix", "groupmix", "groupmix", "uniform", "ints", "offset", "wide", "tiny", "huge", "const", "gappy", "text", "mixed", "empty"]
GOODS = GOOD + ("groupmix",)
LABELSETS = [list("abcdefghijklmnopqrstuvwxy"), ["g%02d" % i for i in range(25)], ["ctl", "trt", "x y", "B", "10", "9", "zz", "é"], ["0", "1", "2", "3", "4", "5", "6"]]

# hand cases of anova: lists of groups of values for column x (None for a missing cell); the group order is the order of the labels as given
FIXED = [
    [[1.0, 2.0, 3.0], [2.0, 3.0, 4.0, 5.0], [7.0, 8.0, 9.5]],
    [[1.5], [2.5, 3.5, 4.0]],
    [[1.0, 2.0], [3.0, 5.0], [4.0, 9.0, 1.0]],
    [[3.0, 3.0, 3.0], [1.0, 2.0, 4.0, 8.0], [5.0, 5.0]],
    [[3.0, 3.0, 3.0], [5.0, 5.0], [7.0, 7.0, 7.0, 7.0]],
    [[3.0, 3.0, 3.0], [3.0, 3.0], [3.0]],
    [[1.0], [2.0], [4.0]],
    [[1.0], [2.0, 3.0, 4.0], [7.0, 7.0, 9.0]],
    [[1e300, 2e300, -3e300, 5e299], [1e300, 3e300, 4e300], [2e300, -1e300, 6e299]],
    [[5e-324, 1e-323, 1.5e-323, 2e-323], [0.0, 5e-324, 1e-323], [3e-323, 1e-323, 4e-323]],
    [[1e8, 100000001.0, 100000002.0, 100000004.0], [100000003.0, 100000007.0, 100000011.0], [100000000.0, 100000006.0, 100000009.0]],
    [[1e8, 100000001.0, 100000002.0], [100000000.5, 100000001.5, 100000002.5], [100000001.0, 100000002.0, 100000003.0]],
    [[1.0, 2.0, 3.0], [], [4.0, 5.0, 7.0]],
    [[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0], [3.0, 4.0, 5.0, 6.0], [1.5, 2.5, 3.5, 4.5]],
    # variances 2^418 apart (Welch: range), the distances of the means 2^637 apart with equal variances (Welch: range), F* near 1e32 (tiny variances, far means), F beyond the largest double (classic)
    [[1.0, 2.0, 3.0], [1e63, 2e63, 3e63], [4.0, 5.0, 6.0]],
    [[-2.0 ** 500, 2.0 ** 500, 0.0], [-2.0 ** 500, 2.0 ** 500, 1e-100], [-2.0 ** 500, 2.0 ** 500, 3.0 * 2.0 ** 302]],
    [[1.0, 1.0000000000000002, 1.0], [2.0, 2.0000000000000004, 2.0], [3.0, 3.0000000000000004, 3.0]],
    [[0.0, 5e-324], [1e300, 1e300], [1e300, 2e300]],
]


def fmt(v):
    return "" if v is None else repr(v)


def run(binary, tmp, names, rows, extra):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "test", "anova", "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=600)


def shuffled(r, seq):
    perm = list(range(len(seq)))
    for i in range(len(perm) - 1, 0, -1):
        k = r.next() % (i + 1)
        perm[i], perm[k] = perm[k], perm[i]
    return [seq[i] for i in perm]


def groupmix(r, labels_of_rows, n):
    """values with a mean and a spread of their own for each group, 3 to 5 digits: the groups differ (or not), as in practice"""
    per = {}
    out = []
    for lab in labels_of_rows:
        if lab not in per:
            per[lab] = (r.uni(-5, 5) * (1 if r.next() % 3 else 0), r.logu(0.05, 20.0))
        m, s = per[lab]
        out.append(repr(round(m + s * (r.u() + r.u() + r.u() - 1.5), 4)) if r.next() % 25 else "")
    return out


class Stat:
    def __init__(self):
        self.ctl = {}
        self.p = {"max": 0.0, "ratio": 0.0, "n": 0}
        self.wd = {"cases": 0, "off": 0, "units": 0}
        self.sp = {"f": 0.0, "p": 0.0, "n": 0, "skipped": 0}
        self.checked = 0
        self.reasons = {}

    def control(self, name, outside):
        c = self.ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0


def ulp_distance(a, b):
    ia = gd.bits(a) if hasattr(gd, "bits") else None
    return None if ia is None else abs(ia - gd.bits(b))


def gen(r, fi):
    """names, rows, group column, the sorted labels, value column indices, groups-of-values per column built later"""
    if fi < 2 * len(FIXED):
        groups = FIXED[fi % len(FIXED)]
        labs = LABELSETS[fi % len(LABELSETS)][:len(groups)]
        rows = []
        for lab, g in zip(labs, groups):
            for v in g:
                rows.append([fmt(v), lab])
            if not g:
                rows.append(["", lab])
        return ["x", "g"], rows, 1, sorted(labs), [0], ["fixed", "group"]
    ls = LABELSETS[int(r.next() % len(LABELSETS))]
    k = [2, 3, 3, 4, 5, 6, 7, 10, 25][int(r.next() % 9)]
    k = min(k, len(ls))
    labels = sorted(ls[:k])
    nval = 1 + int(r.next() % 3)
    shapes = [SHAPES[int(r.next() % len(SHAPES))] for _ in range(nval)]
    n = [k, k + 1, 2 * k, 3 * k, 8, 20, 60, 200][int(r.next() % 8)]
    n = max(n, k)
    grp = []
    weights = [0.2 + r.u() for _ in range(k)]
    tot = sum(weights)
    for i in range(n):
        if r.next() % 100 < 3:
            grp.append("")
            continue
        u = r.u() * tot
        acc = 0.0
        for lab, wgt in zip(labels, weights):
            acc += wgt
            if u <= acc:
                grp.append(lab)
                break
        else:
            grp.append(labels[-1])
    for i in range(min(k, n)):
        grp[i] = labels[i]   # every level has a row (a level with no value in a column comes from the missing cells)
    cols = [groupmix(r, grp, n) if s == "groupmix" else gd.column(r, s, n) for s in shapes]
    gcol = int(r.next() % (nval + 1))
    names = ["c%d_%s" % (j, s) for j, s in enumerate(shapes)]
    allcols = cols[:]
    names.insert(gcol, "g")
    allcols.insert(gcol, grp)
    rows = [[allcols[j][i] for j in range(len(allcols))] for i in range(n)]
    vidx = [j for j in range(len(allcols)) if j != gcol]
    sh = shapes[:]
    sh.insert(gcol, "group")
    return names, rows, gcol, labels, vidx, sh


def column_numeric(rows, j):
    try:
        vals = [parse(rw[j]) for rw in rows]
    except ValueError:
        return False
    return any(v is not None for v in vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scipy", action="store_true")
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    st = Stat()
    files = 0
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            names, rows, gcol, labels, vidx, shapes = gen(r, fi)
            n = len(rows)
            num = [j for j in vidx if column_numeric(rows, j)]
            if not num:
                continue
            chosen = [j for j in num if r.next() % 3 != 0] or [num[0]]
            method = r.pick(["welch", "pooled"]) if fi >= 2 * len(FIXED) else ("welch" if fi < len(FIXED) else "pooled")
            pooled = method == "pooled"
            extra = ["--cols", ",".join(names[j] for j in chosen), "--group", "g", "--variance", method]
            correction = r.pick(["none", "none", "holm", "bh", "bonferroni"])
            if correction != "none":
                extra += ["--correct", correction]
            try:
                p = run(a.bin, tmp, names, rows, extra)
            except subprocess.TimeoutExpired:
                print("gate_anova: the program did not finish (hang)")
                return 2
            if p.returncode != 0:
                keep = pathlib.Path(tempfile.gettempdir()) / ("gate_anova_failed_%d.csv" % fi)
                keep.write_text((tmp / "f.csv").read_text())
                print("gate_anova: exit %d on file %d: %s (the file is %s, args %r)" % (p.returncode, fi, p.stdout[:240], keep, extra))
                return 2
            files += 1
            doc = json.loads(p.stdout)["data"]
            if doc.get("test") != "anova" or doc.get("rows") != n or doc.get("group") != "g" or doc.get("levels") != labels or doc.get("variance") != method or doc.get("adjustment") != correction:
                bad.append((fi, "the header: want anova rows %d group g levels %r variance %s; got %r" % (n, labels, method, {k: doc.get(k) for k in ("test", "rows", "group", "levels", "variance", "adjustment")})))
            tests = doc["tests"]
            if [x["column"] for x in tests] != [names[j] for j in chosen]:
                bad.append((fi, "the tests are %r" % [x["column"] for x in tests]))
                continue
            if correction != "none":
                want_adj = ot.adjust([x["p_value"] for x in tests], correction)
                fam = doc.get("family") or {}
                if fam.get("correction") != correction or fam.get("size") != len([x for x in tests if x["p_value"] is not None]) or not all(gd.same(a_, b_) for a_, b_ in zip(fam.get("p_adjusted", []), want_adj)):
                    bad.append((fi, "the family: want %s %r got %r" % (correction, want_adj, fam)))
            for g, j in zip(tests, chosen):
                gv = [[v for v in (parse(rw[j]) for rw in rows if rw[gcol] == lab) if v is not None] for lab in labels]
                w = ot.anova(gv, not pooled)
                other = ot.anova(gv, pooled)
                if not pooled and w.get("spread", 0) > 368 and w["reason"] is None:
                    # the weights or the distances of the means span more than 2^360: the program says so (a documented limit), the zone 352 to 368 is judged either way
                    st.wd["range"] = st.wd.get("range", 0) + 1
                    if g.get("reason") != "welch-range" or g.get("statistic") is not None or g.get("p_value") is not None:
                        bad.append((fi, "%s: spread %d bits, want welch-range with null statistic, got %r %r" % (tag, w["spread"], g.get("reason"), g.get("statistic"))))
                    continue
                if not pooled and 352 <= w.get("spread", 0) <= 368 and g.get("reason") == "welch-range":
                    st.wd["range"] = st.wd.get("range", 0) + 1
                    continue
                tag = "anova %s %s k=%d" % (names[j], method, len(labels))
                st.checked += 1
                st.reasons[w["reason"]] = st.reasons.get(w["reason"], 0) + 1
                if g["n"] != w["ns"]:
                    bad.append((fi, "%s n: want %r got %r" % (tag, w["ns"], g["n"])))
                for key, idx in (("mean", 1), ("sd", 2)):
                    want = [x[idx] for x in w["groups"]]
                    got = g[key]
                    if len(got) != len(want) or any(not gd.same(x, y) for x, y in zip(got, want)):
                        bad.append((fi, "%s %s: want %r got %r" % (tag, key, want, got)))
                if g.get("method") != method:
                    bad.append((fi, "%s method %r" % (tag, g.get("method"))))
                if g.get("p_method") != "f":
                    bad.append((fi, "%s p_method %r" % (tag, g.get("p_method"))))
                # the statistic: bytes
                if w["statistic"] is None:
                    if g.get("statistic") is not None:
                        bad.append((fi, "%s statistic: want null got %r" % (tag, g.get("statistic"))))
                elif g.get("statistic") is None:
                    bad.append((fi, "%s statistic: want %r got null" % (tag, w["statistic"])))
                elif not gd.same(g["statistic"], w["statistic"]):
                    if not pooled:
                        st.wd["off"] += 1
                        if abs(g["statistic"] - w["statistic"]) > 2 * math.ulp(w["statistic"]) or w["statistic"] == 0.0:
                            bad.append((fi, "%s statistic: want %r got %r" % (tag, w["statistic"], g["statistic"])))
                        else:
                            st.wd["units"] = max(st.wd["units"], round(abs(g["statistic"] - w["statistic"]) / math.ulp(w["statistic"])))
                    else:
                        bad.append((fi, "%s statistic: want %r got %r" % (tag, w["statistic"], g["statistic"])))
                if not pooled and w["statistic"] is not None:
                    st.wd["cases"] += 1
                # degrees of freedom
                df = g.get("df")
                if w["reason"] == "too-few-rows" or (not pooled and w["statistic"] is None):
                    if df is not None:
                        bad.append((fi, "%s df: want null got %r" % (tag, df)))
                elif pooled:
                    if df != {"between": w["k"] - 1, "within": w["N"] - w["k"]}:
                        bad.append((fi, "%s df: want %r got %r" % (tag, {"between": w["k"] - 1, "within": w["N"] - w["k"]}, df)))
                else:
                    if df is None or df["between"] != w["k"] - 1 or not gd.same(df["within"], ot.rounded(w["df2"])):
                        bad.append((fi, "%s df: want between %d within %r got %r" % (tag, w["k"] - 1, ot.rounded(w["df2"]), df)))
                # the effect sizes: bytes
                es = g.get("effect_size")
                if w["eta"] is None:
                    if es is not None:
                        bad.append((fi, "%s effect_size: want null got %r" % (tag, es)))
                elif es is None or es.get("name") != "eta_squared" or not gd.same(es.get("value"), w["eta"]):
                    bad.append((fi, "%s eta: want %r got %r" % (tag, w["eta"], es)))
                if not gd.same(g.get("omega_squared"), w["omega"]):
                    bad.append((fi, "%s omega: want %r got %r" % (tag, w["omega"], g.get("omega_squared"))))
                if g.get("reason") != w["reason"]:
                    bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
                # the p-value
                if w["p"] is None:
                    if g.get("p_value") is not None:
                        bad.append((fi, "%s p_value: want null got %r" % (tag, g.get("p_value"))))
                else:
                    bound = float(w["k"] - 1) + float(w["df2"])
                    if g.get("p_value") is None:
                        bad.append((fi, "%s p_value: want %s got null (reason %r, statistic %r)" % (tag, mpmath.nstr(w["p"], 17), g.get("reason"), g.get("statistic"))))
                        continue
                    if not ot.p_ok(g.get("p_value"), w["p"], bound):
                        bad.append((fi, "%s p_value: want %s got %r" % (tag, mpmath.nstr(w["p"], 17), g.get("p_value"))))
                    elif w["p"] > mpmath.mpf("1e-300"):
                        err = float(abs(mpmath.mpf(g["p_value"]) - w["p"]) / w["p"])
                        st.p["max"] = max(st.p["max"], err)
                        st.p["ratio"] = max(st.p["ratio"], err / (1e-11 + (bound + 1) * ot.EPS))
                        st.p["n"] += 1
                    if not (0.0 <= g["p_value"] <= 1.0):
                        bad.append((fi, "%s p_value %r is not in [0, 1]" % (tag, g["p_value"])))
                    # controls
                    if w["p"] > mpmath.mpf("1e-250") and w["f"] != 0:
                        d1, d2, f = w["k"] - 1, w["df2"], w["f"]
                        sw = ot.p_f(d2 if d2 == int(d2) else d2, d1, f) if d2 == int(d2) else None
                        if sw is not None and d1 != d2:
                            st.control("df-swapped", not ot.p_ok(float(sw), w["p"], bound))
                        st.control("df2-plus-one", not ot.p_ok(float(ot.p_f(d1, d2 + 1, f)), w["p"], bound))
                        if other["p"] is not None:
                            st.control("other-variance-assumption", not ot.p_ok(float(other["p"]), w["p"], bound))
                        chi = mpmath.gammainc(mpmath.mpf(d1) / 2, ot.frac_mpf(d1 * f) / 2, mpmath.inf, regularized=True)
                        st.control("chi-square-approximation", not ot.p_ok(float(chi), w["p"], bound))
                if a.scipy and pooled and w["statistic"] is not None and shapes[j] in GOODS and w["p"] is not None and all(len(x) > 0 for x in gv):
                    import scipy.stats as ss
                    res = ss.f_oneway(*gv)
                    df_ = abs(res.statistic - g["statistic"]) / max(abs(res.statistic), 1.0)
                    dp = abs(res.pvalue - g["p_value"]) / max(abs(res.pvalue), 1e-300)
                    st.sp["n"] += 1
                    st.sp["f"] = max(st.sp["f"], df_)
                    st.sp["p"] = max(st.sp["p"], dp)
                    if df_ > 1e-10 + 1e-14 * abs(res.statistic) or dp > 1e-8:
                        bad.append((fi, "%s scipy: dF %g dp %g" % (tag, df_, dp)))
            if any(b[0] == fi for b in bad):
                pathlib.Path(tempfile.gettempdir(), "gate_anova_bad_%d.csv" % fi).write_text((tmp / "f.csv").read_text() if False else ",".join(names) + "\n" + "\n".join(",".join(rw) for rw in rows) + "\n")
                pathlib.Path(tempfile.gettempdir(), "gate_anova_bad_%d.args" % fi).write_text(" ".join(extra))
            # invariance
            if n > 1:
                p2 = run(a.bin, tmp, names, shuffled(r, rows), extra)
                if p2.stdout != p.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
                if len(chosen) > 1:
                    p3 = run(a.bin, tmp, names, rows, ["--cols", ",".join(reversed([names[j] for j in chosen]))] + extra[2:])
                    if p3.stdout != p.stdout:
                        bad.append((fi, "--cols in another order gives different bytes"))
                # the levels renamed so that the byte order is reversed
                k = len(labels)
                width = len(str(k))
                new = {lab: "q%0*d" % (width, k - 1 - i) for i, lab in enumerate(labels)}
                rrows = [[(new[c] if jj == gcol and c in new else c) for jj, c in enumerate(rw)] for rw in rows]
                p4 = run(a.bin, tmp, names, rrows, extra)
                if p4.returncode != 0:
                    bad.append((fi, "levels renamed: exit %d" % p4.returncode))
                else:
                    for g0, g4 in zip(tests, json.loads(p4.stdout)["data"]["tests"]):
                        ok = g0["p_value"] is None and g4["p_value"] is None or g0["p_value"] is not None and g4["p_value"] is not None
                        ok = ok and g0.get("reason") == g4.get("reason") and g0["n"] == list(reversed(g4["n"])) and g0["mean"] == list(reversed(g4["mean"])) and g0["sd"] == list(reversed(g4["sd"]))
                        ok = ok and g0.get("effect_size") == g4.get("effect_size") and g0.get("omega_squared") == g4.get("omega_squared")
                        if pooled:
                            ok = ok and g0["statistic"] == g4["statistic"] and g0["p_value"] == g4["p_value"] and g0["df"] == g4["df"]
                        elif g0["statistic"] is not None:
                            ok = ok and abs(g0["statistic"] - g4["statistic"]) <= 2 * math.ulp(g0["statistic"]) and abs(g0["p_value"] - g4["p_value"]) <= 1e-9 * max(g0["p_value"], 1e-300)
                        else:
                            ok = ok and g4["statistic"] is None
                        if not ok:
                            bad.append((fi, "the levels renamed do not give the same answer for %s: %r against %r" % (g0["column"], g0, g4)))
    for k, (outside, total) in sorted(st.ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_anova: %d files, %d tests compared (%s), %d differ; p-values against mpmath: largest relative error %.3g over %d, at most %.3g of the bound; Welch statistic: %d compared, %d not the correctly rounded value (at most %d units); controls outside the bound: %s%s" % (
        files, st.checked, ", ".join("%s %d" % (k or "answered", v) for k, v in sorted(st.reasons.items(), key=lambda x: str(x[0]))), len(bad), st.p["max"], st.p["n"], st.p["ratio"], st.wd["cases"], st.wd["off"], st.wd["units"],
        ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(st.ctl.items())),
        ("; scipy f_oneway on %d well-conditioned columns: max relative dF %.3g, dp %.3g" % (st.sp["n"], st.sp["f"], st.sp["p"])) if a.scipy else ""))
    print("  Welch tests refused as welch-range (spread above 2^360): %d" % st.wd.get("range", 0))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
