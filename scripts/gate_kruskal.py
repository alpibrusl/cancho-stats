#!/usr/bin/env python3
"""S4.2 gates R-A, R-B, R-C, R-F and R-G (docs/rank.md section 5) for `stats test kruskal`:

    python3 scripts/gate_kruskal.py [--bin build/stats] [--files 200] [--seed 1] [--scipy]

Generated CSV files (the shapes of gate_describe.py plus lattices with heavy ties, a group column of 2 to 25 levels whose sizes go down to 0, 1 and 2, missing cells in the values and in the group) against
scripts/oracle_test.py, which ranks `Fraction`s by comparing them (not the program's order keys) and shares no code with the program. For each test: n and mean rank of every group (bytes), the number of tie groups,
the statistic H* (tie-corrected: one rounding of an exact rational: bytes), epsilon^2 (bytes), df, the p-value (relative 1e-11 + (k + 1) 2^-53 against mpmath at the EXACT rational statistic), the reasons. Also: the rows
shuffled give the same bytes; the columns named in another order give the same bytes; a monotone transformation of the values (times 8, negated) gives the same statistic, p-value and effect size; the levels renamed so
that their byte order is reversed give the same statistic and the lists reversed. CONTROLS from the same data must leave the bounds on at least a tenth of the cases (of those that can tell): H without the tie
correction (where there are ties), chi-square on k and on k - 2 degrees of freedom, the one-sided p-value. With --scipy: `kruskal` agrees to 1e-12 relative on the statistic and 1e-9 on the p-value on the
well-conditioned columns.

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
from gate_anova import LABELSETS, fmt, shuffled, groupmix, column_numeric  # noqa: E402

SHAPES = ["lattice", "lattice", "lattice", "groupmix", "groupmix", "uniform", "ints", "offset", "wide", "tiny", "huge", "const", "gappy", "text", "mixed", "empty"]
GOODS = ("lattice", "groupmix", "uniform", "ints")

FIXED = [
    [[1.0, 2.0, 3.0, 3.0], [2.0, 3.0, 4.0, 5.0], [7.0, 8.0, 9.5, 3.0]],
    [[1.5], [2.5, 3.5, 4.0]],
    [[1.0], [2.0], [4.0]],
    [[1.0, 2.0], [3.0, 5.0], [4.0, 9.0, 1.0]],
    [[3.0, 3.0, 3.0], [3.0, 3.0], [3.0]],
    [[3.0, 3.0, 3.0], [1.0, 2.0, 4.0, 8.0], [5.0, 5.0]],
    [[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [3.0, 3.0, 3.0]],
    [[1e300, 2e300, -3e300, 5e299], [1e300, 3e300, 4e300], [2e300, -1e300, 6e299]],
    [[5e-324, 1e-323, 1.5e-323, 2e-323], [0.0, 5e-324, 1e-323], [3e-323, 1e-323, 4e-323]],
    [[1.0, 2.0, 3.0], [], [4.0, 5.0, 7.0]],
    [[1.0, 2.0, 3.0], [], []],
    [[1.0, 4.0], [2.0, 3.0]],
    [[1.0, 6.0], [2.0, 5.0], [3.0, 4.0]],
    [[0.0, -0.0, 1.0], [0.0, 2.0], [-1.0, 0.0]],
    [[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0], [3.0, 4.0, 5.0, 6.0], [1.5, 2.5, 3.5, 4.5]],
    [[float(i) for i in range(1, 30)], [float(i) + 0.5 for i in range(1, 30)], [float(i) + 15 for i in range(1, 30)]],
]


def run(binary, tmp, names, rows, extra):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "test", "kruskal", "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=600)


class Stat:
    def __init__(self):
        self.ctl = {}
        self.p = {"max": 0.0, "ratio": 0.0, "n": 0}
        self.sp = {"f": 0.0, "p": 0.0, "n": 0}
        self.checked = 0
        self.reasons = {}

    def control(self, name, outside):
        c = self.ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0


def lattice(r, labels_of_rows, n):
    top = 2 + int(r.next() % 6)
    return ["" if r.next() % 25 == 0 else str(int(r.next() % top)) for _ in range(n)]


def gen(r, fi):
    if fi < len(FIXED):
        groups = FIXED[fi]
        labs = LABELSETS[fi % len(LABELSETS)][:len(groups)]
        rows = []
        for lab, g in zip(labs, groups):
            for v in g:
                rows.append([fmt(v), lab])
            if not g:
                rows.append(["", lab])
        return ["x", "g"], rows, 1, sorted(labs), [0], ["fixed", "group"]
    ls = LABELSETS[int(r.next() % len(LABELSETS))]
    k = min([2, 3, 3, 4, 5, 6, 7, 10, 25][int(r.next() % 9)], len(ls))
    labels = sorted(ls[:k])
    nval = 1 + int(r.next() % 3)
    shapes = [SHAPES[int(r.next() % len(SHAPES))] for _ in range(nval)]
    n = max([k, k + 1, 2 * k, 3 * k, 8, 20, 60, 200][int(r.next() % 8)], k)
    weights = [0.2 + r.u() for _ in range(k)]
    tot = sum(weights)
    grp = []
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
        grp[i] = labels[i]
    cols = []
    for s in shapes:
        if s == "groupmix":
            cols.append(groupmix(r, grp, n))
        elif s == "lattice":
            cols.append(lattice(r, grp, n))
        else:
            cols.append(gd.column(r, s, n))
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
            extra = ["--cols", ",".join(names[j] for j in chosen), "--group", "g"]
            correction = r.pick(["none", "none", "holm", "bh", "bonferroni"])
            if correction != "none":
                extra += ["--correct", correction]
            try:
                p = run(a.bin, tmp, names, rows, extra)
            except subprocess.TimeoutExpired:
                print("gate_kruskal: the program did not finish (hang)")
                return 2
            if p.returncode != 0:
                keep = pathlib.Path(tempfile.gettempdir()) / ("gate_kruskal_failed_%d.csv" % fi)
                keep.write_text((tmp / "f.csv").read_text())
                print("gate_kruskal: exit %d on file %d: %s (the file is %s, args %r)" % (p.returncode, fi, p.stdout[:240], keep, extra))
                return 2
            files += 1
            doc = json.loads(p.stdout)["data"]
            if doc.get("test") != "kruskal" or doc.get("rows") != n or doc.get("group") != "g" or doc.get("levels") != labels or doc.get("adjustment") != correction or "variance" in doc:
                bad.append((fi, "the header: want kruskal rows %d group g levels %r; got %r" % (n, labels, {k: doc.get(k) for k in ("test", "rows", "group", "levels", "variance", "adjustment")})))
            tests = doc["tests"]
            if [x["column"] for x in tests] != [names[j] for j in chosen]:
                bad.append((fi, "the tests are %r" % [x["column"] for x in tests]))
                continue
            if correction != "none":
                want_adj = ot.adjust([x["p_value"] for x in tests], correction)
                fam = doc.get("family") or {}
                if fam.get("correction") != correction or not all(gd.same(a_, b_) for a_, b_ in zip(fam.get("p_adjusted", []), want_adj)):
                    bad.append((fi, "the family: want %s %r got %r" % (correction, want_adj, fam)))
            for g, j in zip(tests, chosen):
                gv = [[v for v in (parse(rw[j]) for rw in rows if rw[gcol] == lab) if v is not None] for lab in labels]
                w = ot.kruskal(gv)
                tag = "kruskal %s k=%d" % (names[j], len(labels))
                st.checked += 1
                st.reasons[w["reason"]] = st.reasons.get(w["reason"], 0) + 1
                if g["n"] != w["ns"]:
                    bad.append((fi, "%s n: want %r got %r" % (tag, w["ns"], g["n"])))
                if len(g["mean_rank"]) != len(w["mean_ranks"]) or any(not gd.same(x, y) for x, y in zip(g["mean_rank"], w["mean_ranks"])):
                    bad.append((fi, "%s mean_rank: want %r got %r" % (tag, w["mean_ranks"], g["mean_rank"])))
                if w["reason"] != "too-few-rows" and g.get("ties") != w["ties"]:
                    bad.append((fi, "%s ties: want %r got %r" % (tag, w["ties"], g.get("ties"))))
                if g.get("method") != "kruskal" or g.get("p_method") != "chi2":
                    bad.append((fi, "%s method %r %r" % (tag, g.get("method"), g.get("p_method"))))
                if not gd.same(g.get("statistic"), w["statistic"]):
                    bad.append((fi, "%s statistic: want %r got %r" % (tag, w["statistic"], g.get("statistic"))))
                if g.get("df") != w["df"]:
                    bad.append((fi, "%s df: want %r got %r" % (tag, w["df"], g.get("df"))))
                es = g.get("effect_size")
                if w["eps"] is None:
                    if es is not None:
                        bad.append((fi, "%s effect_size: want null got %r" % (tag, es)))
                elif es is None or es.get("name") != "epsilon_squared" or not gd.same(es.get("value"), w["eps"]):
                    bad.append((fi, "%s epsilon: want %r got %r" % (tag, w["eps"], es)))
                if g.get("reason") != w["reason"]:
                    bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
                if w["p"] is None:
                    if g.get("p_value") is not None:
                        bad.append((fi, "%s p_value: want null got %r" % (tag, g.get("p_value"))))
                    continue
                bound = float(w["k"])
                if g.get("p_value") is None:
                    bad.append((fi, "%s p_value: want %s got null" % (tag, mpmath.nstr(w["p"], 17))))
                    continue
                if not ot.p_ok(g["p_value"], w["p"], bound):
                    bad.append((fi, "%s p_value: want %s got %r" % (tag, mpmath.nstr(w["p"], 17), g["p_value"])))
                elif w["p"] > mpmath.mpf("1e-300"):
                    err = float(abs(mpmath.mpf(g["p_value"]) - w["p"]) / w["p"])
                    st.p["max"] = max(st.p["max"], err)
                    st.p["ratio"] = max(st.p["ratio"], err / (1e-11 + (bound + 1) * ot.EPS))
                    st.p["n"] += 1
                if w["p"] > mpmath.mpf("1e-250") and w["h"] != 0:
                    d1, h = w["df"], w["h"]
                    nn = w["N"]
                    # H without the tie correction: only where there are ties
                    if w["ties"] > 0:
                        t3 = sum(t ** 3 - t for t in _tie_sizes(gv))
                        h0 = h * (1 - Fraction(t3, nn ** 3 - nn))
                        st.control("no-tie-correction", not ot.p_ok(float(ot.chi2_sf_exact(d1, h0)), w["p"], bound))
                    st.control("chi-square-on-k", not ot.p_ok(float(ot.chi2_sf_exact(d1 + 1, h)), w["p"], bound))
                    if d1 > 1:
                        st.control("chi-square-on-k-minus-2", not ot.p_ok(float(ot.chi2_sf_exact(d1 - 1, h)), w["p"], bound))
                    st.control("one-sided", not ot.p_ok(float(w["p"] / 2), w["p"], bound))
                if a.scipy and shapes[j] in GOODS and all(len(x) > 0 for x in gv) and w["reason"] is None:
                    import scipy.stats as ss
                    res = ss.kruskal(*gv)
                    d_ = abs(res.statistic - g["statistic"]) / max(abs(res.statistic), 1.0)
                    dp = abs(res.pvalue - g["p_value"]) / max(abs(res.pvalue), 1e-300)
                    st.sp["n"] += 1
                    st.sp["f"] = max(st.sp["f"], d_)
                    st.sp["p"] = max(st.sp["p"], dp)
                    if d_ > 1e-12 or dp > 1e-9:
                        bad.append((fi, "%s scipy: dH %g dp %g" % (tag, d_, dp)))
            if any(b[0] == fi for b in bad):
                pathlib.Path(tempfile.gettempdir(), "gate_kruskal_bad_%d.csv" % fi).write_text(",".join(names) + "\n" + "\n".join(",".join(rw) for rw in rows) + "\n")
                pathlib.Path(tempfile.gettempdir(), "gate_kruskal_bad_%d.args" % fi).write_text(" ".join(extra))
            if n > 1:
                p2 = run(a.bin, tmp, names, shuffled(r, rows), extra)
                if p2.stdout != p.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
                if len(chosen) > 1:
                    p3 = run(a.bin, tmp, names, rows, ["--cols", ",".join(reversed([names[j] for j in chosen]))] + extra[2:])
                    if p3.stdout != p.stdout:
                        bad.append((fi, "--cols in another order gives different bytes"))
                # a monotone transformation: times 8 (exact) and negated
                for tf, nm in ((lambda v: v * 8.0, "times 8"), (lambda v: -v, "negated")):
                    trows = []
                    for rw in rows:
                        rw2 = list(rw)
                        for j in chosen:
                            v = parse(rw[j]) if rw[j] != "" else None
                            if v is not None:
                                rw2[j] = repr(tf(v))
                        trows.append(rw2)
                    try:
                        p5 = run(a.bin, tmp, names, trows, extra)
                    except subprocess.TimeoutExpired:
                        bad.append((fi, "transformed file: hang"))
                        continue
                    if p5.returncode != 0:
                        # an overflow of the transformed values is not the program's fault: the huge shapes times 8
                        if not any(shapes[j] in ("huge", "wide") for j in chosen):
                            bad.append((fi, "transformed (%s) file: exit %d" % (nm, p5.returncode)))
                        continue
                    for g0, g5 in zip(tests, json.loads(p5.stdout)["data"]["tests"]):
                        keys = ("statistic", "p_value", "effect_size", "df", "ties", "n", "reason")
                        if any(g0.get(k_) != g5.get(k_) for k_ in keys):
                            if not any(shapes[j] in ("huge",) for j in chosen):
                                bad.append((fi, "a monotone transformation (%s) changes the answer: %r against %r" % (nm, {k_: g0.get(k_) for k_ in keys}, {k_: g5.get(k_) for k_ in keys})))
                k = len(labels)
                width = len(str(k))
                new = {lab: "q%0*d" % (width, k - 1 - i) for i, lab in enumerate(labels)}
                rrows = [[(new[c] if jj == gcol and c in new else c) for jj, c in enumerate(rw)] for rw in rows]
                p4 = run(a.bin, tmp, names, rrows, extra)
                if p4.returncode != 0:
                    bad.append((fi, "levels renamed: exit %d" % p4.returncode))
                else:
                    for g0, g4 in zip(tests, json.loads(p4.stdout)["data"]["tests"]):
                        ok = all(g0.get(k_) == g4.get(k_) for k_ in ("statistic", "p_value", "effect_size", "df", "ties", "reason"))
                        ok = ok and g0["n"] == list(reversed(g4["n"])) and g0["mean_rank"] == list(reversed(g4["mean_rank"]))
                        if not ok:
                            bad.append((fi, "the levels renamed do not give the same answer for %s: %r against %r" % (g0["column"], g0, g4)))
    for k, (outside, total) in sorted(st.ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_kruskal: %d files, %d tests compared (%s), %d differ; p-values against mpmath: largest relative error %.3g over %d, at most %.3g of the bound; controls outside the bound: %s%s" % (
        files, st.checked, ", ".join("%s %d" % (k or "answered", v) for k, v in sorted(st.reasons.items(), key=lambda x: str(x[0]))), len(bad), st.p["max"], st.p["n"], st.p["ratio"],
        ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(st.ctl.items())),
        ("; scipy kruskal on %d well-conditioned columns: max relative dH %.3g, dp %.3g" % (st.sp["n"], st.sp["f"], st.sp["p"])) if a.scipy else ""))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


def _tie_sizes(gv):
    allv = [Fraction(v) for g in gv for v in g]
    out = []
    for v in set(allv):
        t = allv.count(v)
        if t > 1:
            out.append(t)
    return out


if __name__ == "__main__":
    sys.exit(main())
