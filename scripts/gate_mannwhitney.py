#!/usr/bin/env python3
"""S4.3 gates R-A, R-B, R-C, R-F and R-G (docs/rank.md section 5) for `stats test mannwhitney`:

    python3 scripts/gate_mannwhitney.py [--bin build/stats] [--files 200] [--seed 1] [--scipy]

Generated CSV files (the shapes of gate_describe.py plus lattices with heavy ties, a group column of exactly two levels whose sizes go down to 0, 1 and 2, missing cells in the values and in the group, sizes up to 200 so that
both sides of the exact limit of 60 are met) against scripts/oracle_test.py, which takes U from the pairs (no ranks), counts the arrangements by an independent recurrence (checked first against the enumeration of every
choice of ranks for every size up to 14) and shares no code with the program. For each test and each of `--exact auto|yes|no` and `--continuity auto|none`: n of each group, the ties, U of both groups, the rank-biserial
correlation, z (bytes), the p-value (the EXACT one equal to the oracle's rational rounded once, bit for bit; the normal one against mpmath at the exact z^2, relative 1e-11 + 2^-53), `p_method`, the reasons; `--exact yes`
where the exact p-value is not available is refused with stats.exact-unavailable. Also: the rows shuffled, `--cols` in another order, the values times 8 give the same bytes; the values negated give U_b for U_a, z and the
rank-biserial correlation negated, the same p; the levels renamed so that their byte order flips do the same. CONTROLS from the same data must leave the bounds on at least a tenth of the cases that can tell: the normal p
without the continuity correction, without the tie correction (where there are ties), the normal p where the exact one is the answer, the one-sided p. With --scipy: `mannwhitneyu` (method exact and asymptotic) agrees to 1e-12
relative on U and on the exact p, 1e-9 on the asymptotic p, on the well-conditioned columns.

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
from gate_test import parse  # noqa: E402
from gate_anova import fmt, shuffled, groupmix, column_numeric  # noqa: E402
from gate_kruskal import lattice  # noqa: E402

SHAPES = ["lattice", "lattice", "groupmix", "groupmix", "groupmix", "uniform", "ints", "offset", "wide", "tiny", "huge", "const", "gappy", "text", "mixed", "empty"]
GOODS = ("lattice", "groupmix", "uniform", "ints")
LABELS = [("a", "b"), ("ctl", "trt"), ("0", "1"), ("10", "9"), ("B", "a"), ("x y", "x z")]

FIXED = [
    ([1.1, 2.3, 3.3, 4.0, 5.5], [2.0, 3.0, 6.0, 7.0, 8.0, 9.0]),
    ([1.0, 2.0, 2.0, 3.0, 5.0], [2.0, 3.0, 3.0, 6.0, 7.0, 8.0]),
    ([1.0], [2.0]),
    ([2.0], [1.0]),
    ([1.0, 2.0], [3.0, 4.0]),
    ([3.0, 3.0, 3.0], [3.0, 3.0]),
    ([1.0, 4.0], [2.0, 3.0]),
    ([5e-324, 1e-323, 1.5e-323], [0.0, 5e-324, 2e-323]),
    ([1e300, -1e300, 2e300], [1e299, 3e300]),
    ([float(i) for i in range(1, 31)], [float(i) + 0.5 for i in range(1, 31)]),
    ([float(i) for i in range(1, 32)], [float(i) + 0.5 for i in range(1, 31)]),
    ([float(i) for i in range(1, 31)], [float(i) + 0.5 for i in range(1, 32)]),
    ([0.0, -0.0, 1.0], [0.0, 2.0]),
    ([1.0, 2.0, 3.0], []),
    ([float(i) for i in range(1, 31)], [float(i) + 0.25 for i in range(1, 31)]),
]


def run(binary, tmp, names, rows, extra):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "test", "mannwhitney", "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=600)


class Stat:
    def __init__(self):
        self.ctl = {}
        self.p = {"max": 0.0, "n": 0}
        self.exact = {"n": 0, "differs": 0}
        self.sp = {"u": 0.0, "pe": 0.0, "pn": 0.0, "n": 0}
        self.checked = 0
        self.refused = 0
        self.reasons = {}

    def control(self, name, outside):
        c = self.ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0


def gen(r, fi):
    if fi < len(FIXED):
        ga, gb = FIXED[fi]
        labels = LABELS[fi % len(LABELS)]
        rows = [[fmt(v), labels[0]] for v in ga] + [[fmt(v), labels[1]] for v in gb]
        if not ga:
            rows.append(["", labels[0]])
        if not gb:
            rows.append(["", labels[1]])
        return ["x", "g"], rows, 1, sorted(labels), [0], ["fixed", "group"]
    labels = sorted(LABELS[int(r.next() % len(LABELS))])
    nval = 1 + int(r.next() % 3)
    shapes = [SHAPES[int(r.next() % len(SHAPES))] for _ in range(nval)]
    n = [2, 3, 4, 5, 8, 12, 20, 40, 60, 61, 62, 200][int(r.next() % 12)]
    pa = [0.5, 0.5, 0.2, 0.8, 0.08][int(r.next() % 5)]
    grp = []
    for i in range(n):
        u = r.u()
        if u < 0.03:
            grp.append("")
        else:
            grp.append(labels[0] if r.u() < pa else labels[1])
    if n >= 2:
        grp[0], grp[1] = labels[0], labels[1]
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
    # the oracle's own counts against the enumeration of every choice of ranks
    for na in range(1, 8):
        for nb in range(1, 8):
            if na + nb <= 14 and ot.mw_counts(na, nb) != ot.mw_brute(na, nb):
                print("gate_mannwhitney: the oracle's counts differ from the enumeration at %d, %d" % (na, nb))
                return 2
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
            exact = r.pick(["auto", "auto", "no", "yes"])
            cont = r.pick(["auto", "auto", "none"])
            extra = ["--cols", ",".join(names[j] for j in chosen), "--group", "g", "--exact", exact, "--continuity", cont]
            correction = r.pick(["none", "none", "holm", "bh"])
            if correction != "none":
                extra += ["--correct", correction]
            gvs = {}
            for j in chosen:
                gvs[j] = [[v for v in (parse(rw[j]) for rw in rows if rw[gcol] == lab) if v is not None] for lab in labels]
            ora = {j: ot.mannwhitney(gvs[j][0], gvs[j][1], exact, cont == "auto") for j in chosen}
            # the group column must have both levels in the file; the generator makes sure of it for n >= 2
            if len(set(rw[gcol] for rw in rows if rw[gcol] != "")) != 2:
                continue
            try:
                p = run(a.bin, tmp, names, rows, extra)
            except subprocess.TimeoutExpired:
                print("gate_mannwhitney: the program did not finish (hang)")
                return 2
            want_refusal = exact == "yes" and any(not ora[j]["exact_available"] and ora[j]["u_a"] is not None for j in chosen)
            # a refusal also when a column has no value in one group? no: that is too-few-rows, a null; the exact check is on the columns with both groups
            if want_refusal:
                st.refused += 1
                doc = json.loads(p.stdout) if p.stdout else {}
                if p.returncode != 2 or doc.get("error", {}).get("rule") != "stats.exact-unavailable":
                    bad.append((fi, "--exact yes where the exact p-value is not available: want exit 2 stats.exact-unavailable, got exit %d %r" % (p.returncode, p.stdout[:160])))
                continue
            if p.returncode != 0:
                keep = pathlib.Path(tempfile.gettempdir()) / ("gate_mannwhitney_failed_%d.csv" % fi)
                keep.write_text((tmp / "f.csv").read_text())
                print("gate_mannwhitney: exit %d on file %d: %s (the file is %s, args %r)" % (p.returncode, fi, p.stdout[:240], keep, extra))
                return 2
            files += 1
            doc = json.loads(p.stdout)["data"]
            if doc.get("test") != "mannwhitney" or doc.get("rows") != n or doc.get("group") != "g" or doc.get("levels") != labels or doc.get("adjustment") != correction:
                bad.append((fi, "the header: want mannwhitney rows %d levels %r; got %r" % (n, labels, {k: doc.get(k) for k in ("test", "rows", "group", "levels", "adjustment")})))
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
                w = ora[j]
                tag = "mannwhitney %s %s %s" % (names[j], exact, cont)
                st.checked += 1
                st.reasons[w["reason"]] = st.reasons.get(w["reason"], 0) + 1
                if g["n"] != {"a": w["na"], "b": w["nb"]}:
                    bad.append((fi, "%s n: want %r got %r" % (tag, {"a": w["na"], "b": w["nb"]}, g["n"])))
                if g.get("reason") != w["reason"]:
                    bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
                if g.get("method") != "mannwhitney" or g.get("continuity") != (cont == "auto"):
                    bad.append((fi, "%s method/continuity %r %r" % (tag, g.get("method"), g.get("continuity"))))
                if w["u_a"] is None:
                    if g.get("statistic") is not None or g.get("p_value") is not None:
                        bad.append((fi, "%s: want nulls got %r" % (tag, g)))
                    continue
                if g["ties"] != w["ties"]:
                    bad.append((fi, "%s ties: want %r got %r" % (tag, w["ties"], g["ties"])))
                nanb = w["na"] * w["nb"]
                ub = Fraction(nanb) - w["u_a"]
                if not gd.same(g.get("statistic"), w["statistic"]) or g["u"] != {"a": w["statistic"], "b": ot.rounded(ub)} and not (gd.same(g["u"]["a"], w["statistic"]) and gd.same(g["u"]["b"], ot.rounded(ub))):
                    bad.append((fi, "%s U: want %r %r got %r %r" % (tag, w["statistic"], ot.rounded(ub), g.get("statistic"), g.get("u"))))
                es = g.get("effect_size")
                if es is None or es.get("name") != "rank_biserial" or not gd.same(es.get("value"), w["rb"]):
                    bad.append((fi, "%s rank_biserial: want %r got %r" % (tag, w["rb"], es)))
                if w["z"] is None:
                    if g.get("z") is not None:
                        bad.append((fi, "%s z: want null got %r" % (tag, g.get("z"))))
                elif not gd.same(g.get("z"), w["z"]):
                    bad.append((fi, "%s z: want %r got %r" % (tag, w["z"], g.get("z"))))
                if w["p"] is None:
                    if g.get("p_value") is not None:
                        bad.append((fi, "%s p_value: want null got %r" % (tag, g.get("p_value"))))
                    continue
                if g.get("p_method") != w["p_method"]:
                    bad.append((fi, "%s p_method: want %r got %r" % (tag, w["p_method"], g.get("p_method"))))
                if g.get("p_value") is None:
                    bad.append((fi, "%s p_value: want %s got null" % (tag, mpmath.nstr(w["p"], 17))))
                    continue
                if w["p_method"] == "exact":
                    st.exact["n"] += 1
                    if not gd.same(g["p_value"], ot.rounded(w["p_exact"])):
                        st.exact["differs"] += 1
                        bad.append((fi, "%s exact p: want %r got %r" % (tag, ot.rounded(w["p_exact"]), g["p_value"])))
                else:
                    if not ot.p_ok(g["p_value"], w["p"], 1.0):
                        bad.append((fi, "%s normal p: want %s got %r" % (tag, mpmath.nstr(w["p"], 17), g["p_value"])))
                    elif w["p"] > mpmath.mpf("1e-300"):
                        st.p["max"] = max(st.p["max"], float(abs(mpmath.mpf(g["p_value"]) - w["p"]) / w["p"]))
                        st.p["n"] += 1
                # controls
                if w["z"] is not None and w["p"] > mpmath.mpf("1e-250"):
                    z2 = w["z2"]
                    if w["p_method"] == "normal" and cont == "auto":
                        nn = w["na"] + w["nb"]
                        base = w["z2"]
                        # without the continuity correction: m' = |m| instead of |m| - 1
                        mm = abs(int(2 * w["u_a"]) - nanb)
                        if mm > 1:
                            t3 = _t3(gvs[j])
                            z0 = Fraction(3 * mm * mm * nn * (nn - 1), nanb * (nn ** 3 - nn - t3))
                            st.control("no-continuity-correction", not ot.p_ok(float(min(mpmath.mpf(1), ot.norm_sf_exact(z0))), w["p"], 1.0))
                    if w["p_method"] == "normal" and w["ties"] > 0:
                        nn = w["na"] + w["nb"]
                        mm = abs(int(2 * w["u_a"]) - nanb) - (1 if cont == "auto" and int(2 * w["u_a"]) != nanb else 0)
                        z0 = Fraction(3 * mm * mm * nn * (nn - 1), nanb * (nn ** 3 - nn))
                        st.control("no-tie-correction", not ot.p_ok(float(min(mpmath.mpf(1), ot.norm_sf_exact(z0))), w["p"], 1.0))
                    if w["p_method"] == "exact" and w["z2"] is not None:
                        st.control("normal-for-exact", not ot.p_ok(float(min(mpmath.mpf(1), ot.norm_sf_exact(w["z2"]))), w["p"], 1.0))
                    st.control("one-sided", not ot.p_ok(float(w["p"] / 2), w["p"], 1.0))
                if a.scipy and shapes[j] in GOODS and w["na"] > 0 and w["nb"] > 0 and w["reason"] is None:
                    import scipy.stats as ss
                    res = ss.mannwhitneyu(gvs[j][0], gvs[j][1], method="exact" if w["p_method"] == "exact" else "asymptotic", use_continuity=(cont == "auto"))
                    du = abs(res.statistic - g["statistic"]) / max(abs(res.statistic), 1.0)
                    dp = abs(res.pvalue - g["p_value"]) / max(abs(res.pvalue), 1e-300)
                    st.sp["n"] += 1
                    st.sp["u"] = max(st.sp["u"], du)
                    key = "pe" if w["p_method"] == "exact" else "pn"
                    st.sp[key] = max(st.sp[key], dp)
                    if du > 1e-12 or dp > (1e-12 if key == "pe" else 1e-9):
                        bad.append((fi, "%s scipy: dU %g dp %g" % (tag, du, dp)))
            if any(b[0] == fi for b in bad):
                pathlib.Path(tempfile.gettempdir(), "gate_mannwhitney_bad_%d.csv" % fi).write_text(",".join(names) + "\n" + "\n".join(",".join(rw) for rw in rows) + "\n")
                pathlib.Path(tempfile.gettempdir(), "gate_mannwhitney_bad_%d.args" % fi).write_text(" ".join(extra))
            if n > 1:
                p2 = run(a.bin, tmp, names, shuffled(r, rows), extra)
                if p2.stdout != p.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
                if len(chosen) > 1:
                    p3 = run(a.bin, tmp, names, rows, ["--cols", ",".join(reversed([names[j] for j in chosen]))] + extra[2:])
                    if p3.stdout != p.stdout:
                        bad.append((fi, "--cols in another order gives different bytes"))
                # the values times 8 (same bytes of every number), negated (the groups swap their order in the ranking)
                for nm, tf in (("times 8", lambda v: v * 8.0), ("negated", lambda v: -v)):
                    trows = []
                    for rw in rows:
                        rw2 = list(rw)
                        for j in chosen:
                            if rw[j] != "":
                                rw2[j] = repr(tf(parse(rw[j])))
                        trows.append(rw2)
                    p5 = run(a.bin, tmp, names, trows, extra)
                    if p5.returncode != 0:
                        if not any(shapes[j] in ("huge", "wide") for j in chosen):
                            bad.append((fi, "transformed (%s) file: exit %d" % (nm, p5.returncode)))
                        continue
                    for g0, g5 in zip(tests, json.loads(p5.stdout)["data"]["tests"]):
                        if g0.get("statistic") is None:
                            continue
                        if nm == "times 8":
                            keys = ("statistic", "u", "z", "p_value", "p_method", "effect_size", "n", "ties", "reason")
                            ok = all(g0.get(k_) == g5.get(k_) for k_ in keys)
                        else:
                            neg = lambda x: None if x is None else (-x if x != 0 else 0.0)
                            ok = g0["p_value"] == g5["p_value"] and g0["p_method"] == g5["p_method"] and g0["n"] == g5["n"] and g0["ties"] == g5["ties"]
                            ok = ok and g0["u"]["a"] == g5["u"]["b"] and g0["u"]["b"] == g5["u"]["a"] and g0["statistic"] == g5["u"]["a"] - 0 if False else ok
                            ok = ok and g0["u"]["a"] == g5["u"]["b"] and g0["u"]["b"] == g5["u"]["a"]
                            ok = ok and gd.same(neg(g0["z"]), g5["z"]) and gd.same(neg(g0["effect_size"]["value"]), g5["effect_size"]["value"])
                        if not ok and not any(shapes[j] in ("huge",) for j in chosen):
                            bad.append((fi, "the values %s change the answer: %r against %r" % (nm, g0, g5)))
                k_ = labels
                new = {labels[0]: "zz" + labels[1], labels[1]: "yy" + labels[0]}
                srows = [[(new[c] if jj == gcol and c in new else c) for jj, c in enumerate(rw)] for rw in rows]
                p4 = run(a.bin, tmp, names, srows, extra)
                if p4.returncode != 0:
                    bad.append((fi, "levels renamed: exit %d" % p4.returncode))
                else:
                    neg = lambda x: None if x is None else (-x if x != 0 else 0.0)
                    for g0, g4 in zip(tests, json.loads(p4.stdout)["data"]["tests"]):
                        if g0.get("statistic") is None:
                            continue
                        ok = g0["p_value"] == g4["p_value"] and g0["p_method"] == g4["p_method"] and g0["ties"] == g4["ties"] and g0["n"] == {"a": g4["n"]["b"], "b": g4["n"]["a"]}
                        ok = ok and g0["u"]["a"] == g4["u"]["b"] and g0["u"]["b"] == g4["u"]["a"] and gd.same(neg(g0["z"]), g4["z"]) and gd.same(neg(g0["effect_size"]["value"]), g4["effect_size"]["value"])
                        if not ok:
                            bad.append((fi, "the groups swapped do not mirror the answer: %r against %r" % (g0, g4)))
    for k, (outside, total) in sorted(st.ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_mannwhitney: %d files, %d tests compared (%s), %d refusals expected and seen, %d differ; exact p-values equal to the oracle's rational rounded once: %d of %d; normal p-values against mpmath: largest relative error %.3g over %d; controls outside the bound: %s%s" % (
        files, st.checked, ", ".join("%s %d" % (k or "answered", v) for k, v in sorted(st.reasons.items(), key=lambda x: str(x[0]))), st.refused, len(bad), st.exact["n"] - st.exact["differs"], st.exact["n"], st.p["max"], st.p["n"],
        ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(st.ctl.items())),
        ("; scipy on %d columns: max relative dU %.3g, exact dp %.3g, normal dp %.3g" % (st.sp["n"], st.sp["u"], st.sp["pe"], st.sp["pn"])) if a.scipy else ""))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


def _t3(groups):
    allv = [Fraction(v) for g in groups for v in g]
    t3 = 0
    for v in set(allv):
        t = allv.count(v)
        if t > 1:
            t3 += t ** 3 - t
    return t3


if __name__ == "__main__":
    sys.exit(main())
