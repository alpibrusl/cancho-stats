#!/usr/bin/env python3
"""S4.5 gates K-A, K-B, K-C and K-F (docs/rank.md section 7) for `stats test ks`:

    python3 scripts/gate_ks.py [--bin build/stats] [--files 200] [--seed 1] [--scipy]

Generated CSV files (the shapes of gate_describe.py plus lattices with heavy ties, a group column of exactly two levels whose sizes go down to 0, 1 and 2, missing cells in the values and in the group, sizes up to 200 so that both
sides of the exact limit of 60 are met) against scripts/oracle_test.py, which takes D from the two empirical distribution functions at every distinct value (Fractions, no sorting of the pooled sample), counts the arrangements
whose distance is at least the observed one by a recurrence over the lattice that counts the paths that have LEFT the band (the program counts those that stay inside; the two are checked against the enumeration of every choice for
every size up to 5 x 5 first), and takes the asymptotic p-value from the alternating series at 160 digits at the exact lambda. For each test and each of `--exact auto|yes|no`: n of each group, the tie groups, d and D (bytes), lambda
(bytes), the p-value (the EXACT one equal to the oracle's rational rounded once, bit for bit; the asymptotic one against the oracle, relative 1e-11 + 2^-53), `p_method`, the reasons; `--exact yes` where the exact p-value is not
available (any two values equal, or more than 60 values) is refused with stats.exact-unavailable. Also: the rows shuffled, `--cols` in another order, the values times 8, the values negated and the levels renamed so that their
byte order flips give the same bytes of D, the p-value and the method (the sizes swapped for the renamed levels). CONTROLS from the same data must leave the bounds on at least a tenth of the cases that can tell: the asymptotic p
where the exact one is the answer, Stephens' corrected lambda, half the p-value. With --scipy: `ks_2samp` (method exact) agrees on D to 1e-12 and on the exact p-value to 1e-12 relative, on the well-conditioned columns (scipy's
asymptotic p-value is a different approximation and is not compared; docs/rank.md section 13).

Exit 0: every number is the expected one. Exit 1: a difference. Exit 2: the harness failed.
"""
import argparse
import json
import math
import pathlib
import subprocess
import sys
import tempfile

import mpmath

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
import gate_describe as gd  # noqa: E402
import oracle_test as ot  # noqa: E402
from gate_test import parse  # noqa: E402
from gate_anova import fmt, shuffled, groupmix, column_numeric  # noqa: E402
from gate_kruskal import lattice  # noqa: E402
from gate_mannwhitney import SHAPES, GOODS, LABELS, FIXED, gen  # noqa: E402


def run(binary, tmp, names, rows, extra):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "test", "ks", "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=600)


def pooled_ok(groups):
    v = [x for g in groups for x in g]
    return len(v) <= 60 and len(set(v)) == len(v)


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
    checked = refused = files = 0
    exact_n = exact_differs = 0
    pmax = 0.0
    pn = 0
    sp = {"d": 0.0, "p": 0.0, "n": 0}
    reasons = {}
    for na in range(1, 6):
        for nb in range(1, 6):
            for d in range(1, na * nb + 1):
                if ot.ks_brute(na, nb, d) != ot.ks_exceed(na, nb, d):
                    print("gate_ks: the oracle's counts differ from the enumeration at %d, %d, %d" % (na, nb, d))
                    return 2

    def control(name, outside):
        c = ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0

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
            extra = ["--cols", ",".join(names[j] for j in chosen), "--group", "g", "--exact", exact]
            correction = r.pick(["none", "none", "holm", "bh"])
            if correction != "none":
                extra += ["--correct", correction]
            gvs = {}
            for j in chosen:
                gvs[j] = [[v for v in (parse(rw[j]) for rw in rows if rw[gcol] == lab) if v is not None] for lab in labels]
            ora = {j: ot.ks(gvs[j][0], gvs[j][1], exact) for j in chosen}
            if len(set(rw[gcol] for rw in rows if rw[gcol] != "")) != 2:
                continue
            try:
                p = run(a.bin, tmp, names, rows, extra)
            except subprocess.TimeoutExpired:
                print("gate_ks: the program did not finish (hang)")
                return 2
            # the refusal is a property of the pooled values of the column (any two equal, or more than 60), whether or not a group is empty
            want_refusal = exact == "yes" and any(not pooled_ok(gvs[j]) for j in chosen)
            if want_refusal:
                refused += 1
                doc = json.loads(p.stdout) if p.stdout else {}
                if p.returncode != 2 or doc.get("error", {}).get("rule") != "stats.exact-unavailable":
                    bad.append((fi, "--exact yes where the exact p-value is not available: want exit 2 stats.exact-unavailable, got exit %d %r" % (p.returncode, p.stdout[:160])))
                continue
            if p.returncode != 0:
                keep = pathlib.Path(tempfile.gettempdir()) / ("gate_ks_failed_%d.csv" % fi)
                keep.write_text((tmp / "f.csv").read_text())
                print("gate_ks: exit %d on file %d: %s (the file is %s, args %r)" % (p.returncode, fi, p.stdout[:240], keep, extra))
                return 2
            files += 1
            doc = json.loads(p.stdout)["data"]
            if doc.get("test") != "ks" or doc.get("rows") != n or doc.get("group") != "g" or doc.get("levels") != labels or doc.get("adjustment") != correction:
                bad.append((fi, "the header: want ks rows %d levels %r; got %r" % (n, labels, {k: doc.get(k) for k in ("test", "rows", "group", "levels", "adjustment")})))
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
                tag = "ks %s %s" % (names[j], exact)
                checked += 1
                reasons[w["reason"]] = reasons.get(w["reason"], 0) + 1
                if g["n"] != {"a": w["na"], "b": w["nb"]}:
                    bad.append((fi, "%s n: want %r got %r" % (tag, {"a": w["na"], "b": w["nb"]}, g["n"])))
                if g.get("reason") != w["reason"]:
                    bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
                if g.get("method") != "ks":
                    bad.append((fi, "%s method %r" % (tag, g.get("method"))))
                if w["d"] is None:
                    if g.get("statistic") is not None or g.get("p_value") is not None or g.get("distance") is not None:
                        bad.append((fi, "%s: want nulls got %r" % (tag, g)))
                    continue
                if g["ties"] != w["ties"]:
                    bad.append((fi, "%s ties: want %r got %r" % (tag, w["ties"], g["ties"])))
                if g.get("distance") != w["d"] or not gd.same(g.get("statistic"), w["statistic"]) or not gd.same(g.get("lambda"), w["lam"]):
                    bad.append((fi, "%s d, D, lambda: want %r %r %r got %r %r %r" % (tag, w["d"], w["statistic"], w["lam"], g.get("distance"), g.get("statistic"), g.get("lambda"))))
                if g.get("p_method") != w["p_method"]:
                    bad.append((fi, "%s p_method: want %r got %r" % (tag, w["p_method"], g.get("p_method"))))
                if g.get("p_value") is None:
                    bad.append((fi, "%s p_value: want %s got null" % (tag, mpmath.nstr(w["p"], 17))))
                    continue
                if w["p_method"] == "exact":
                    exact_n += 1
                    if not gd.same(g["p_value"], ot.rounded(w["p_exact"])):
                        exact_differs += 1
                        bad.append((fi, "%s exact p: want %r got %r" % (tag, ot.rounded(w["p_exact"]), g["p_value"])))
                else:
                    if not ot.p_ok(g["p_value"], w["p"], 1.0):
                        bad.append((fi, "%s asymptotic p: want %s got %r" % (tag, mpmath.nstr(w["p"], 17), g["p_value"])))
                    elif w["p"] > mpmath.mpf("1e-300"):
                        pmax = max(pmax, float(abs(mpmath.mpf(g["p_value"]) - w["p"]) / w["p"]))
                        pn += 1
                # controls
                if w["d"] > 0 and w["p"] > mpmath.mpf("1e-250"):
                    nn = w["na"] + w["nb"]
                    if w["p_method"] == "exact":
                        x = mpmath.sqrt(ot.frac_mpf(w["lam2"]))
                        q = sum(2 * (-1) ** (k - 1) * mpmath.exp(-2 * k * k * x * x) for k in range(1, 200))
                        control("asymptotic-for-exact", not ot.p_ok(float(min(mpmath.mpf(1), q)), w["p"], 1.0))
                    else:
                        # Stephens: lambda (sqrt(en) + 0.12 + 0.11 / sqrt(en)) / sqrt(en)
                        en = mpmath.mpf(w["na"] * w["nb"]) / nn
                        x = mpmath.sqrt(ot.frac_mpf(w["lam2"])) * (mpmath.sqrt(en) + mpmath.mpf("0.12") + mpmath.mpf("0.11") / mpmath.sqrt(en)) / mpmath.sqrt(en)
                        q = sum(2 * (-1) ** (k - 1) * mpmath.exp(-2 * k * k * x * x) for k in range(1, 200))
                        control("stephens", not ot.p_ok(float(min(mpmath.mpf(1), q)), w["p"], 1.0))
                    control("half-p", not ot.p_ok(float(w["p"] / 2), w["p"], 1.0))
                if a.scipy and shapes[j] in GOODS and w["na"] > 0 and w["nb"] > 0 and w["p_method"] == "exact":
                    import scipy.stats as ss
                    res = ss.ks_2samp(gvs[j][0], gvs[j][1], method="exact")
                    dd = abs(res.statistic - g["statistic"]) / max(abs(res.statistic), 1e-300)
                    dp = abs(res.pvalue - g["p_value"]) / max(abs(res.pvalue), 1e-300)
                    sp["n"] += 1
                    sp["d"] = max(sp["d"], dd)
                    sp["p"] = max(sp["p"], dp)
                    if dd > 1e-12 or dp > 1e-12:
                        bad.append((fi, "%s scipy: dD %g dp %g" % (tag, dd, dp)))
            if any(b[0] == fi for b in bad):
                pathlib.Path(tempfile.gettempdir(), "gate_ks_bad_%d.csv" % fi).write_text(",".join(names) + "\n" + "\n".join(",".join(rw) for rw in rows) + "\n")
                pathlib.Path(tempfile.gettempdir(), "gate_ks_bad_%d.args" % fi).write_text(" ".join(extra))
            if n > 1:
                p2 = run(a.bin, tmp, names, shuffled(r, rows), extra)
                if p2.stdout != p.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
                if len(chosen) > 1:
                    p3 = run(a.bin, tmp, names, rows, ["--cols", ",".join(reversed([names[j] for j in chosen]))] + extra[2:])
                    if p3.stdout != p.stdout:
                        bad.append((fi, "--cols in another order gives different bytes"))
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
                        keys = ("statistic", "distance", "lambda", "p_value", "p_method", "n", "ties", "reason")
                        if not all(g0.get(k_) == g5.get(k_) for k_ in keys) and not any(shapes[j] in ("huge",) for j in chosen):
                            bad.append((fi, "the values %s change the answer: %r against %r" % (nm, g0, g5)))
                new = {labels[0]: "zz" + labels[1], labels[1]: "yy" + labels[0]}
                srows = [[(new[c] if jj == gcol and c in new else c) for jj, c in enumerate(rw)] for rw in rows]
                p4 = run(a.bin, tmp, names, srows, extra)
                if p4.returncode != 0:
                    bad.append((fi, "levels renamed: exit %d" % p4.returncode))
                else:
                    for g0, g4 in zip(tests, json.loads(p4.stdout)["data"]["tests"]):
                        ok = all(g0.get(k_) == g4.get(k_) for k_ in ("statistic", "distance", "lambda", "p_value", "p_method", "ties", "reason")) and g0["n"] == {"a": g4["n"]["b"], "b": g4["n"]["a"]}
                        if not ok:
                            bad.append((fi, "the groups swapped do not mirror the answer: %r against %r" % (g0, g4)))
    for k, (outside, total) in sorted(ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_ks: %d files, %d tests compared (%s), %d refusals expected and seen, %d differ; exact p-values equal to the oracle's rational rounded once: %d of %d; asymptotic p-values against mpmath: largest relative error %.3g over %d; controls outside the bound: %s%s" % (
        files, checked, ", ".join("%s %d" % (k or "answered", v) for k, v in sorted(reasons.items(), key=lambda x: str(x[0]))), refused, len(bad), exact_n - exact_differs, exact_n, pmax, pn,
        ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(ctl.items())),
        ("; scipy on %d columns: max relative dD %.3g, exact dp %.3g" % (sp["n"], sp["d"], sp["p"])) if a.scipy else ""))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
