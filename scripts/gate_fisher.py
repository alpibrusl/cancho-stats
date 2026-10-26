#!/usr/bin/env python3
"""S4.6 gates F-A, F-B, F-C and F-F (docs/rank.md section 14) for `stats test fisher`:

    python3 scripts/gate_fisher.py [--bin build/stats] [--files 200] [--seed 1] [--scipy]

Generated CSV files (two to five categorical columns of two levels with seeded skew and dependence, labels whose byte order differs from their first appearance, missing cells, sizes from 2 to 3,000, and tables built from exact
counts: the empty cells, the margins of 0, N = 10,000 and N = 10,001, and the TOLERANCE TABLES found by search, where another table has a probability larger than the observed one by a relative 1.1e-8 to 1.5e-8, so that
the tolerance of 1e-7 in R and scipy would add it to the sum and the definition does not) against scripts/oracle_test.py, which sums `math.comb` weights as Fractions (checked first against the enumeration of every choice
of the first column for every table up to 12 items) and shares no code with the program. For each test: the table, n and missing, the odds ratio (bytes) or its reason (`odds-ratio-undefined`, `odds-ratio-infinite`), and the
p-value EQUAL to the oracle's rational rounded once, bit for bit; a column with a level count other than 2 (in it or in --vs) is refused with stats.group-levels, a table of more than 10,000 rows with limit.exact-too-large
(exit 8). Also: the rows shuffled, `--cols` in another order, the labels renamed so that their byte order flips (the odds ratio inverted, the p-value the same), the two columns exchanged (the same odds ratio and p-value) give
the same numbers. CONTROLS from the same data must leave the bound on at least a tenth of the cases that can tell: the p-value with R's tolerance of 1e-7 (judged on the tolerance tables), twice the smaller tail, the
one-sided tail. With --scipy: `fisher_exact` agrees on the odds ratio (1e-12 relative) and on the p-value (1e-12 relative) wherever the tolerance does not matter (the oracle's two sums are equal).

Exit 0: every number is the expected one. Exit 1: a difference. Exit 2: the harness failed.
"""
import argparse
import json
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
import gate_describe as gd  # noqa: E402
import oracle_test as ot  # noqa: E402
from gate_anova import shuffled  # noqa: E402

LABELS = [("a", "b"), ("no", "yes"), ("0", "1"), ("10", "9"), ("B", "a"), ("x y", "x z"), ("ctl", "trt")]

# (a, b, c, d): the first row is the first level of the column, the first column the first level of --vs
FIXED = [
    (3, 1, 1, 4), (0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 1, 1), (0, 5, 5, 0), (5, 0, 0, 5), (10, 0, 0, 0), (0, 3, 0, 4), (2, 0, 3, 0), (0, 2, 3, 0), (4, 0, 0, 0), (1, 2, 3, 4),
    (50, 79, 159, 202),     # N = 490: the table with a = 60 has a probability larger by 1.1e-8 relative
    (60, 69, 149, 212),
    (131, 601, 169, 598),   # N = 1499: a = 162 larger by 1.5e-8
    (162, 570, 138, 629),
    (12, 316, 96, 639),     # N = 1063: a = 57 larger by 4.9e-8
    (57, 271, 51, 684),
    (5000, 0, 0, 5000), (2500, 2500, 2500, 2500), (3000, 2000, 1000, 4000), (9999, 0, 0, 1), (4000, 1, 1, 5998), (10000, 1, 0, 0), (5000, 1, 4999, 5), (1, 0, 2, 1),
]
TOLERANCE = {(50, 79, 159, 202), (131, 601, 169, 598), (12, 316, 96, 639)}


def table_rows(r, tab, lx, ly):
    a, b, c, d = tab
    rows = []
    for (i, j), k in (((0, 0), a), ((0, 1), b), ((1, 0), c), ((1, 1), d)):
        rows += [[lx[i], ly[j]] for _ in range(k)]
    return rows


def run(binary, tmp, names, rows, extra, kind="fisher"):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "test", kind, "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=600)


def counts(rows, j, k):
    """the table of the complete rows of the columns j and k: the levels by byte order, the 2 x 2 counts or None when a column does not have two levels, and the number of complete rows"""
    comp = [(rw[j], rw[k]) for rw in rows if rw[j] != "" and rw[k] != ""]
    lj = sorted(set(x for x, _ in comp))
    lk = sorted(set(y for _, y in comp))
    return lj, lk, comp


def gen(r, fi):
    if fi < len(FIXED):
        lx = sorted(LABELS[fi % len(LABELS)])
        ly = sorted(LABELS[(fi * 3 + 1) % len(LABELS)])
        rows = table_rows(r, FIXED[fi], lx, ly)
        return ["x", "y"], rows, [0]
    ncol = 2 + int(r.next() % 4)
    n = [2, 3, 5, 8, 12, 20, 40, 100, 300, 1000, 3000][int(r.next() % 11)]
    labs = [sorted(LABELS[int(r.next() % len(LABELS))]) for _ in range(ncol)]
    skew = [0.5, 0.5, 0.2, 0.8, 0.05, 0.95][int(r.next() % 6)]
    dep = [0.0, 0.3, 0.7][int(r.next() % 3)]
    rows = []
    for _ in range(n):
        y = 1 if r.u() < skew else 0
        row = []
        for j in range(ncol - 1):
            x = y if r.u() < dep else (1 if r.u() < 0.5 else 0)
            row.append("" if r.u() < 0.03 else labs[j][x])
        row.append("" if r.u() < 0.03 else labs[ncol - 1][y])
        rows.append(row)
    if n >= 2:
        for j in range(ncol):
            rows[0][j] = labs[j][0]
            rows[1][j] = labs[j][1]
    names = ["c%d" % j for j in range(ncol - 1)] + ["y"]
    return names, rows, list(range(ncol - 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scipy", action="store_true")
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    ctl = {}
    checked = refused = files = exact_n = exact_differs = 0
    sp = {"o": 0.0, "p": 0.0, "n": 0}

    def control(name, outside):
        c = ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0

    import itertools
    for t in itertools.product(range(4), repeat=4):
        n = sum(t)
        if 0 < n <= 12:
            a_, b_, c_, d_ = t
            r1, r2, c1 = a_ + b_, c_ + d_, a_ + c_
            import math
            if ot.fisher_brute(*t) != {x: math.comb(r1, x) * math.comb(r2, c1 - x) for x in range(max(0, c1 - r2), min(r1, c1) + 1)}:
                print("gate_fisher: the oracle's weights differ from the enumeration at %r" % (t,))
                return 2
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            names, rows, vidx = gen(r, fi)
            n = len(rows)
            chosen = [j for j in vidx if fi < len(FIXED) or r.next() % 3 != 0] or [vidx[0]]
            extra = ["--cols", ",".join(names[j] for j in chosen), "--vs", "y"]
            correction = r.pick(["none", "none", "holm", "bh"])
            if correction != "none":
                extra += ["--correct", correction]
            ycol = len(names) - 1
            # the expected outcome, column by column in the order of --cols
            expected = []
            refusal = None
            for j in chosen:
                lj, lk, comp = counts(rows, j, ycol)
                if not comp:
                    expected.append((j, None, lj, lk))
                    continue
                if len(lj) != 2 or len(lk) != 2:
                    refusal = ("stats.group-levels", 2)
                    break
                if len(comp) > 10000:
                    refusal = ("limit.exact-too-large", 8)
                    break
                cnt = {(i, k): 0 for i in range(2) for k in range(2)}
                for x, y in comp:
                    cnt[(lj.index(x), lk.index(y))] += 1
                expected.append((j, ot.fisher(cnt[(0, 0)], cnt[(0, 1)], cnt[(1, 0)], cnt[(1, 1)]), lj, lk))
            try:
                p = run(a.bin, tmp, names, rows, extra)
            except subprocess.TimeoutExpired:
                print("gate_fisher: the program did not finish (hang)")
                return 2
            if refusal is not None:
                refused += 1
                doc = json.loads(p.stdout) if p.stdout else {}
                if p.returncode != refusal[1] or doc.get("error", {}).get("rule") != refusal[0]:
                    bad.append((fi, "want exit %d %s, got exit %d %r" % (refusal[1], refusal[0], p.returncode, p.stdout[:160])))
                continue
            if p.returncode != 0:
                keep = pathlib.Path(tempfile.gettempdir()) / ("gate_fisher_failed_%d.csv" % fi)
                keep.write_text((tmp / "f.csv").read_text())
                print("gate_fisher: exit %d on file %d: %s (the file is %s, args %r)" % (p.returncode, fi, p.stdout[:240], keep, extra))
                return 2
            files += 1
            doc = json.loads(p.stdout)["data"]
            if doc.get("test") != "fisher" or doc.get("vs") != "y" or doc.get("adjustment") != correction:
                bad.append((fi, "the header %r" % {k: doc.get(k) for k in ("test", "vs", "adjustment")}))
            tests = doc["tests"]
            if [x["column"] for x in tests] != [names[j] for j in chosen]:
                bad.append((fi, "the tests are %r" % [x["column"] for x in tests]))
                continue
            if correction != "none":
                want_adj = ot.adjust([x["p_value"] for x in tests], correction)
                fam = doc.get("family") or {}
                if fam.get("correction") != correction or not all(gd.same(a_, b_) for a_, b_ in zip(fam.get("p_adjusted", []), want_adj)):
                    bad.append((fi, "the family: want %s %r got %r" % (correction, want_adj, fam)))
            for g, (j, w, lj, lk) in zip(tests, expected):
                tag = "fisher %s" % names[j]
                checked += 1
                if w is None or w["status"] is not None:
                    if g.get("reason") != "too-few-rows" or g.get("p_value") is not None or g.get("statistic") is not None:
                        bad.append((fi, "%s: want too-few-rows got %r" % (tag, g)))
                    continue
                if g["n"]["n"] != w["n"] or g["n"]["missing"] != n - w["n"]:
                    bad.append((fi, "%s n: want %d missing %d got %r" % (tag, w["n"], n - w["n"], g["n"])))
                if g.get("method") != "fisher" or g.get("interval") is not None or g.get("df") is not None:
                    bad.append((fi, "%s method/interval/df %r" % (tag, {k: g.get(k) for k in ("method", "interval", "df")})))
                if g["assumptions"]["levels"] != {"column": lj, "vs": lk}:
                    bad.append((fi, "%s levels %r" % (tag, g["assumptions"]["levels"])))
                if w["odds"] is None:
                    if g.get("statistic") is not None or g.get("effect_size") is not None or g.get("reason") != w["odds_reason"]:
                        bad.append((fi, "%s: want odds null %s got %r" % (tag, w["odds_reason"], g)))
                else:
                    want = ot.rounded(w["odds"])
                    es = g.get("effect_size") or {}
                    if not gd.same(g.get("statistic"), want) or es.get("name") != "odds_ratio" or not gd.same(es.get("value"), want) or g.get("reason") is not None:
                        bad.append((fi, "%s odds ratio: want %r got %r %r reason %r" % (tag, want, g.get("statistic"), es, g.get("reason"))))
                wp = ot.rounded(w["p"])
                exact_n += 1
                if not gd.same(g.get("p_value"), wp) or g.get("p_method") != "exact, probabilities not larger than the observed":
                    exact_differs += 1
                    bad.append((fi, "%s p: want %r got %r (%r)" % (tag, wp, g.get("p_value"), g.get("p_method"))))
                # controls
                if (g_tab := g["assumptions"]["table"]) and tuple(g_tab) in TOLERANCE:
                    control("tolerance-1e-7", not gd.same(g.get("p_value"), ot.rounded(w["p_tol"])))
                if w["n"] >= 4:
                    control("twice-the-smaller-tail", not gd.same(g.get("p_value"), ot.rounded(w["p_double"])))
                    control("one-sided-tail", not gd.same(g.get("p_value"), ot.rounded(w["p_less"])))
                if a.scipy and w["odds"] is not None and w["p"] == w["p_tol"]:
                    import scipy.stats as ss
                    tab = list(g["assumptions"]["table"])
                    res = ss.fisher_exact([[tab[0], tab[1]], [tab[2], tab[3]]])
                    do = abs(res.statistic - g["statistic"]) / max(abs(res.statistic), 1e-300)
                    dp = abs(res.pvalue - g["p_value"]) / max(abs(res.pvalue), 1e-300)
                    sp["n"] += 1
                    sp["o"] = max(sp["o"], do)
                    sp["p"] = max(sp["p"], dp)
                    if do > 1e-12 or dp > 1e-12:
                        bad.append((fi, "%s scipy: dodds %g dp %g" % (tag, do, dp)))
            if any(b[0] == fi for b in bad):
                pathlib.Path(tempfile.gettempdir(), "gate_fisher_bad_%d.csv" % fi).write_text(",".join(names) + "\n" + "\n".join(",".join(rw) for rw in rows) + "\n")
                pathlib.Path(tempfile.gettempdir(), "gate_fisher_bad_%d.args" % fi).write_text(" ".join(extra))
            if n > 1 and n <= 3000:
                p2 = run(a.bin, tmp, names, shuffled(r, rows), extra)
                if p2.stdout != p.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
                if len(chosen) > 1:
                    p3 = run(a.bin, tmp, names, rows, ["--cols", ",".join(reversed([names[j] for j in chosen]))] + extra[2:])
                    if correction == "none":
                        # the tests follow the order of --cols; each test is the same
                        d3 = {x["column"]: x for x in json.loads(p3.stdout)["data"]["tests"]}
                        if p3.returncode != 0 or any(d3.get(x["column"]) != x for x in tests):
                            bad.append((fi, "--cols in another order changes a test"))
                # the labels of the column renamed so that the byte order flips: the rows swap, the odds ratio is inverted, the p-value stays
                j0 = chosen[0]
                lj0 = sorted(set(rw[j0] for rw in rows if rw[j0] != ""))
                if len(lj0) == 2 and correction == "none":
                    new = {lj0[0]: "zz" + lj0[1], lj0[1]: "yy" + lj0[0]}
                    srows = [[(new[c] if jj == j0 and c in new else c) for jj, c in enumerate(rw)] for rw in rows]
                    p4 = run(a.bin, tmp, names, srows, ["--cols", names[j0], "--vs", "y"])
                    p5 = run(a.bin, tmp, names, rows, ["--cols", names[j0], "--vs", "y"])
                    if p4.returncode == 0 and p5.returncode == 0:
                        g4 = json.loads(p4.stdout)["data"]["tests"][0]
                        g5 = json.loads(p5.stdout)["data"]["tests"][0]
                        if g5.get("p_value") != g4.get("p_value") or g5["assumptions"].get("table") and g4["assumptions"]["table"] != [g5["assumptions"]["table"][2], g5["assumptions"]["table"][3], g5["assumptions"]["table"][0], g5["assumptions"]["table"][1]]:
                            bad.append((fi, "the rows swapped do not mirror the table or keep the p-value: %r against %r" % (g5, g4)))
                    # the two columns exchanged
                    p6 = run(a.bin, tmp, names, rows, ["--cols", "y", "--vs", names[j0]])
                    if p6.returncode == 0 and p5.returncode == 0:
                        g5 = json.loads(p5.stdout)["data"]["tests"][0]
                        g6 = json.loads(p6.stdout)["data"]["tests"][0]
                        if g5.get("p_value") != g6.get("p_value") or g5.get("statistic") != g6.get("statistic"):
                            bad.append((fi, "exchanging the columns changes the odds ratio or the p-value: %r against %r" % (g5, g6)))
    for k, (outside, total) in sorted(ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_fisher: %d files, %d tests compared, %d refusals expected and seen, %d differ; p-values equal to the oracle's rational rounded once: %d of %d; controls outside the bound: %s%s" % (
        files, checked, refused, len(bad), exact_n - exact_differs, exact_n, ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(ctl.items())),
        ("; scipy on %d tables: max relative dodds %.3g, dp %.3g" % (sp["n"], sp["o"], sp["p"])) if a.scipy else ""))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
