#!/usr/bin/env python3
"""S3.3 gates T-A, T-B, T-D, T-H and T-I (docs/test.md sections 5 and 10) for `stats test chisq`:

    python3 scripts/gate_test3.py [--bin build/stats] [--files 200] [--seed 1] [--scipy]

Generated CSV files with categorical columns (labels that are text, numbers, with spaces, in any order of appearance; empty cells in every column; tables with empty cells, one level, tall and wide) against
scripts/oracle_test.py, which works from `Fraction` and `mpmath` and shares no code with the program. Independence (`--cols A,B --vs G`) and goodness of fit (`--cols A,B --expected ..`, proportions of 1 to 9 decimals
drawn to add up to exactly 1): the statistic against the CORRECTLY ROUNDED exact value (bytes; the program sums in double-double, so a mismatch is counted and reported: it is expected to be 0 and to be possible only
within 2^-100 of a tie), `df`, `n` and the missing rows, the numbers of expected counts below 5 (exactly) and the smallest expected count (2 ulp), the effect size (4 ulp of the square root of the exact value), the
p-value (relative 1e-11 + (X2 + 1) 2^-53 against mpmath at the exact X2), the table (or the levels in the order of their bytes), the reasons. Also: the rows shuffled give the same bytes; the columns named in another
order give the same entries in that order. CONTROLS must leave the bound on a tenth of the cases: Yates's correction (2 x 2 and larger, applied to every cell), the degrees of freedom r c - 1 (k for a fit), the G
statistic in place of Pearson's. With --scipy: `chi2_contingency(correction=False)` and `chisquare` agree to 1e-12 (statistic) and 1e-9 (p-value) relative.

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

POOLS = [["a", "b", "c", "d", "e", "f", "g"], ["1", "2", "3", "10", "9", "100", "20"], ["red", "green", "blue", "x y", "Z", "yes", "no"], ["0", "1", "1.0", "01", "-1", "2e0", "NA"], ["L%02d" % i for i in range(40)], ["x,1", "", "y", "z w", "q\"r", "-"]]
ULP = 2.0 ** -52


def raw(v):
    """the bytes the program sees between the quotes of a quoted cell: a doubled quote stays doubled (as `describe --by` names its groups), so a level is named by its raw text and ordered by it"""
    return v.replace('"', '""')


def csv_cell(v):
    """None is an empty cell (missing); the empty label is written as a quoted empty cell, which is a level; a label with a comma or a quote is quoted (RFC 4180)"""
    if v is None:
        return ""
    if v == "" or "," in v or '"' in v:
        return '"' + v.replace('"', '""') + '"'
    return v


def run(binary, tmp, names, rows, extra):
    text = ",".join(names) + "\n" + "\n".join(",".join(csv_cell(c) for c in r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "test", "chisq", "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=300)


def draw_table(r):
    """(labels of the first column, labels of the second, the weights of the cells)"""
    pa = POOLS[int(r.next() % len(POOLS))]
    pb = POOLS[int(r.next() % len(POOLS))]
    ra = [1, 2, 2, 3, 3, 4, 5, 7, 7, 7, 20, 40][int(r.next() % 12)]
    cb = [1, 2, 2, 2, 3, 3, 4, 6, 7, 7, 15, 40][int(r.next() % 12)]
    la = pa[:ra]
    lb = pb[:cb]
    w = [[0 if r.next() % 7 == 0 else 1 + int(r.next() % 20) for _ in lb] for _ in la]
    return la, lb, w


def sample(r, la, lb, w, n):
    cells = [(i, j) for i in range(len(la)) for j in range(len(lb))]
    tot = sum(w[i][j] for i, j in cells)
    out = []
    for _ in range(n):
        if tot == 0:
            out.append((la[0], lb[0]))
            continue
        t = int(r.next() % tot)
        for i, j in cells:
            t -= w[i][j]
            if t < 0:
                out.append((la[i], lb[j]))
                break
    return out


def draw_props(r, k):
    """k proportions with 1 to 9 decimals that add up to exactly 1, each above 0 (as the integers over 10^d and the decimal texts)"""
    d = 1 + int(r.next() % 9)
    scale = 10 ** d
    if k > scale:
        d = 3
        scale = 1000
    wts = [1 + int(r.next() % 20) for _ in range(k)]
    tot = sum(wts)
    parts = [max(1, scale * wi // tot) for wi in wts]
    parts[0] += scale - sum(parts)
    if min(parts) < 1:
        parts = [scale // k] * k
        parts[0] += scale - sum(parts)
    texts = []
    for p in parts:
        s = "%0*d" % (d + 1, p)
        texts.append(s[:-d] + "." + s[-d:])
    return [Fraction(p, scale) for p in parts], texts


def ulps_ok(got, want, k):
    if got is None or want is None:
        return got is None and want is None
    if want == 0:
        return got == 0
    return abs(got - want) <= k * ULP * abs(want)


def g_pvalue(df, obs, exp):
    g = mpmath.mpf(0)
    for o, e in zip(obs, exp):
        if o > 0:
            g += 2 * o * mpmath.log(mpmath.mpf(o) / ot.frac_mpf(e))
    return mpmath.gammainc(mpmath.mpf(df) / 2, g / 2, mpmath.inf, regularized=True)


class Stat:
    def __init__(self):
        self.ctl = {}
        self.p = {"max": 0.0, "ratio": 0.0, "n": 0}
        self.sp = {"x": 0.0, "p": 0.0, "n": 0}
        self.checked = 0
        self.tie = 0
        self.effect = 0.0

    def control(self, name, outside):
        c = self.ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0


def check_entry(bad, fi, tag, g, w, fit, rows_total, levels, st):
    st.checked += 1
    if g.get("method") != ("goodness of fit" if fit else "independence") or g.get("rounding") != "double-double" or g.get("p_method") != "chi2" or g.get("interval") is not None:
        bad.append((fi, "%s method/rounding/p_method/interval: %r" % (tag, {k: g.get(k) for k in ("method", "rounding", "p_method", "interval")})))
    if g["n"] != {"n": w["n"], "missing": rows_total - w["n"]}:
        bad.append((fi, "%s n: want %r got %r" % (tag, {"n": w["n"], "missing": rows_total - w["n"]}, g["n"])))
    if g.get("reason") != w["reason"]:
        bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
    if w["statistic"] is None:
        for k in ("statistic", "df", "p_value", "effect_size"):
            if g.get(k) is not None:
                bad.append((fi, "%s %s: want null got %r" % (tag, k, g.get(k))))
    else:
        if not gd.same(g.get("statistic"), w["statistic"]):
            st.tie += 1
            bad.append((fi, "%s statistic: want %r (the correctly rounded exact value) got %r" % (tag, w["statistic"], g.get("statistic"))))
        if g.get("df") != w["df"] or not isinstance(g.get("df"), int):
            bad.append((fi, "%s df: want %d got %r" % (tag, w["df"], g.get("df"))))
        es = g.get("effect_size")
        name = "cohen_w" if fit else "cramers_v"
        if es is None or es.get("name") != name or es.get("rounding") != "double" or not ulps_ok(es.get("value"), float(w["effect"]), 4):
            bad.append((fi, "%s effect_size: want %s %s got %r" % (tag, name, mpmath.nstr(w["effect"], 17), es)))
        else:
            st.effect = max(st.effect, float(abs(mpmath.mpf(es["value"]) - w["effect"]) / (ULP * w["effect"]))) if w["effect"] != 0 else st.effect
        if not ot.p_ok(g.get("p_value"), w["p"], float(w["x2"])):
            bad.append((fi, "%s p_value: want %s got %r" % (tag, mpmath.nstr(w["p"], 17), g.get("p_value"))))
        elif w["p"] > mpmath.mpf("1e-300"):
            err = float(abs(mpmath.mpf(g["p_value"]) - w["p"]) / w["p"])
            st.p["max"] = max(st.p["max"], err)
            st.p["ratio"] = max(st.p["ratio"], err / (1e-11 + (float(w["x2"]) + 1) * ot.EPS))
            st.p["n"] += 1
    a = g.get("assumptions") or {}
    if w["n"] > 0:
        if a.get("expected_below_5") != w["below5"]:
            bad.append((fi, "%s expected_below_5: want %d got %r" % (tag, w["below5"], a.get("expected_below_5"))))
        if not ulps_ok(a.get("min_expected"), float(w["min_expected"]), 2):
            bad.append((fi, "%s min_expected: want %r got %r" % (tag, float(w["min_expected"]), a.get("min_expected"))))
        if fit:
            if a.get("levels") != levels:
                bad.append((fi, "%s levels: want %r got %r" % (tag, levels, a.get("levels"))))
        elif a.get("table") != {"rows": w["rows"], "cols": w["cols"]}:
            bad.append((fi, "%s table: want %r got %r" % (tag, {"rows": w["rows"], "cols": w["cols"]}, a.get("table"))))
    elif a:
        bad.append((fi, "%s assumptions of an empty table: %r" % (tag, a)))


def controls(st, g, w, fit):
    if w["statistic"] is None or w["p"] is None or not w["p"] > mpmath.mpf("1e-250") or w["df"] < 1:
        return
    x2 = float(w["x2"])
    df = w["df"]
    obs = w["obs"]
    exp = w["expected_flat"]
    if not fit:
        wrong = ot.chi2_sf(w["rows"] * w["cols"] - 1, w["x2"])
        st.control("df-rc-minus-one", not ot.p_ok(float(wrong), w["p"], x2))
        y = sum((abs(o - e) - Fraction(1, 2)) ** 2 / e for o, e in zip(obs, exp) if abs(o - e) >= Fraction(1, 2))
        y = sum(((abs(o - e) - Fraction(1, 2)) ** 2 if abs(o - e) > Fraction(1, 2) else Fraction(0)) / e for o, e in zip(obs, exp))
        st.control("yates", not ot.p_ok(float(ot.chi2_sf(df, y)), w["p"], x2))
    else:
        wrong = ot.chi2_sf(len(obs), w["x2"])
        st.control("df-k", not ot.p_ok(float(wrong), w["p"], x2))
    st.control("g-statistic", not ot.p_ok(float(g_pvalue(df, obs, exp)), w["p"], x2))


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
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            la, lb, wt = draw_table(r)
            n = [2, 3, 5, 10, 30, 100, 400, 1500][int(r.next() % 8)]
            if fi % 50 == 49:
                n = 120000   # N times a count above 2^31: the integers of the double-double terms are large
            pairs = sample(r, la, lb, wt, n)
            # a second and a third categorical column, drawn independently of the first
            lc = POOLS[int(r.next() % len(POOLS))][:1 + int(r.next() % 5)]
            third = [lc[int(r.next() % len(lc))] for _ in range(n)]
            cols = {"A": [p[0] for p in pairs], "G": [p[1] for p in pairs], "C": third}
            # empty cells
            for name in cols:
                for i in range(n):
                    if r.next() % 20 == 0:
                        cols[name][i] = None
            names = ["A", "G", "C"]
            rows = [[cols[nm][i] for nm in names] for i in range(n)]
            # --- independence of A and C against G
            extra = ["--cols", "A,C", "--vs", "G"]
            p = run(a.bin, tmp, names, rows, extra)
            if p.returncode != 0:
                keep = pathlib.Path(tempfile.gettempdir()) / ("gate_test3_failed_%d.csv" % fi)
                keep.write_text((tmp / "f.csv").read_text())
                print("gate_test3: exit %d on file %d: %s (the file is %s, args %r)" % (p.returncode, fi, p.stdout[:240], keep, extra))
                return 2
            doc = json.loads(p.stdout)["data"]
            if doc.get("test") != "chisq" or doc.get("vs") != "G" or doc.get("adjustment") != "none" or [x["column"] for x in doc["tests"]] != ["A", "C"]:
                bad.append((fi, "the header: %r" % {k: doc.get(k) for k in ("test", "vs", "adjustment")}))
                continue
            for g, nm in zip(doc["tests"], ["A", "C"]):
                used = [(cols[nm][i], cols["G"][i]) for i in range(n) if cols[nm][i] is not None and cols["G"][i] is not None]
                ra = sorted(set(x for x, _ in used), key=raw)
                cb = sorted(set(y for _, y in used), key=raw)
                table = [[sum(1 for u in used if u == (x, y)) for y in cb] for x in ra]
                w = ot.chisq_independence(table) if used else ot.chisq_independence([])
                if used:
                    w["obs"] = [v for row in table for v in row]
                    w["expected_flat"] = [e for row in w["expected"] for e in row]
                check_entry(bad, fi, "chisq %s vs G" % nm, g, w, False, n, None, st)
                controls(st, g, w, False)
                if a.scipy and w["statistic"] is not None and g.get("p_value") is not None and used:
                    import scipy.stats as ss
                    res = ss.chi2_contingency(table, correction=False)
                    dx = abs(res[0] - g["statistic"]) / max(abs(res[0]), 1e-300)
                    dp = abs(res[1] - g["p_value"]) / max(abs(res[1]), 1e-300)
                    st.sp["n"] += 1
                    st.sp["x"] = max(st.sp["x"], dx)
                    st.sp["p"] = max(st.sp["p"], dp)
                    if dx > 1e-12 or dp > 1e-9:
                        bad.append((fi, "scipy chi2_contingency: dx %g dp %g" % (dx, dp)))
            # --- the rows shuffled, the columns in another order
            if n > 1:
                perm = list(range(n))
                for i in range(n - 1, 0, -1):
                    k = r.next() % (i + 1)
                    perm[i], perm[k] = perm[k], perm[i]
                p2 = run(a.bin, tmp, names, [rows[i] for i in perm], extra)
                if p2.stdout != p.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
                p3 = run(a.bin, tmp, names, rows, ["--cols", "C,A", "--vs", "G"])
                d3 = json.loads(p3.stdout)["data"]["tests"] if p3.returncode == 0 else None
                if d3 is None or d3 != list(reversed(doc["tests"])):
                    bad.append((fi, "--cols in another order does not give the same entries in that order"))
            # --- goodness of fit of A and C
            used_levels = {}
            for nm in ("A", "C"):
                vals = [v for v in cols[nm] if v is not None]
                used_levels[nm] = sorted(set(vals), key=raw)
            if len(set(len(used_levels[nm]) for nm in used_levels)) == 1 and used_levels["A"]:
                k = len(used_levels["A"])
                props, texts = draw_props(r, k)
                extra = ["--cols", "A,C", "--expected", ",".join(texts)]
                p = run(a.bin, tmp, names, rows, extra)
                if p.returncode != 0:
                    print("gate_test3: exit %d on a fit of file %d: %s (args %r)" % (p.returncode, fi, p.stdout[:240], extra))
                    return 2
                doc = json.loads(p.stdout)["data"]
                if doc.get("test") != "chisq" or doc.get("expected") != [float(x) for x in props]:
                    bad.append((fi, "the header of the fit: %r" % {k2: doc.get(k2) for k2 in ("test", "expected")}))
                for g, nm in zip(doc["tests"], ["A", "C"]):
                    vals = [v for v in cols[nm] if v is not None]
                    lv = used_levels[nm]
                    counts = [sum(1 for v in vals if v == x) for x in lv]
                    w = ot.chisq_fit(counts, props)
                    w["obs"] = counts
                    w["expected_flat"] = w["expected"]
                    check_entry(bad, fi, "chisq %s fit" % nm, g, w, True, n, [raw(x) for x in lv], st)
                    controls(st, g, w, True)
                    if a.scipy and w["statistic"] is not None and g.get("p_value") is not None:
                        import scipy.stats as ss
                        res = ss.chisquare(counts, [float(e) for e in w["expected"]])
                        dx = abs(res.statistic - g["statistic"]) / max(abs(res.statistic), 1e-300)
                        dp = abs(res.pvalue - g["p_value"]) / max(abs(res.pvalue), 1e-300)
                        st.sp["n"] += 1
                        st.sp["x"] = max(st.sp["x"], dx)
                        st.sp["p"] = max(st.sp["p"], dp)
                        if dx > 1e-12 or dp > 1e-9:
                            bad.append((fi, "scipy chisquare: dx %g dp %g" % (dx, dp)))
    for k, (outside, total) in sorted(st.ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_test3 chisq: %d files, %d entries compared, %d differ (statistics not equal to the correctly rounded exact value: %d); p-values against mpmath: largest relative error %.3g over %d, at most %.3g of the bound; effect size within %.3g ulp of the exact root; controls outside the bound: %s%s" % (
        a.files, st.checked, len(bad), st.tie, st.p["max"], st.p["n"], st.p["ratio"], st.effect, ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(st.ctl.items())),
        ("; scipy on %d tables: max relative dX2 %.3g, dp %.3g" % (st.sp["n"], st.sp["x"], st.sp["p"])) if a.scipy else ""))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
