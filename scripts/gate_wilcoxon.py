#!/usr/bin/env python3
"""S4.4 gates R-A, R-B, R-C, R-F and R-G (docs/rank.md section 5) for `stats test wilcoxon`:

    python3 scripts/gate_wilcoxon.py [--bin build/stats] [--files 200] [--seed 1] [--scipy]

Generated CSV files (pairs of columns that share a part, lattices with zero differences and heavy ties, offsets, wide exponents, missing cells; the one-sample form against --mu0 too; sizes up to 200 so that both sides of the exact
limit of 50 nonzero differences are met) against scripts/oracle_test.py, which takes the differences as EXACT Fractions, ranks |d| by comparing them (no order keys, no two-sum) and shares no code with the program; the exact counts are
by the product of (1 + q^i), checked against the enumeration of every subset up to n = 8 (14 in the unit test of the oracle). For each test, each of the three methods for zeros (`wilcox`, `pratt`, `zsplit`), `--exact auto|yes|no`,
`--continuity auto|none`: the counts of rows, nonzero differences and zeros, the ties, W+ and W-, the rank-biserial correlation, z (bytes), the p-value (the EXACT one equal to the oracle's rational rounded once, bit for bit; the
normal one against mpmath at the exact z^2, relative 1e-11 + 2^-53), `p_method`, the reasons; `--exact yes` where the exact p-value is not available is refused (stats.exact-unavailable). Also: the rows shuffled and `--cols` /
`--pairs` reordered give the same bytes; the values times 8 give the same bytes; the pair named the other way round (and the values negated) exchange W+ and W-, negate z and the rank-biserial correlation, keep the p-value.
CONTROLS from the same data must leave the bounds on at least a tenth of the cases that can tell: the normal p without the continuity correction, without the tie correction (where there are ties), the normal p where the exact
one is the answer, the p of another method for zeros (where there are zeros), the one-sided p. With --scipy: `wilcoxon` agrees on the statistic (min(W+, W-)) and the p-values (exact 1e-12, normal 1e-9 relative) for the three
methods for zeros where every difference is exact in floating point (scipy subtracts in floating point: a difference that rounds can create a tie that the exact differences do not have, and the program ranks the exact ones).

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
from gate_anova import fmt, shuffled, column_numeric  # noqa: E402

GOODS = ("lattice", "shared", "uniform", "ints")
MU0S = ["0", "1", "0.5", "2.5", "-1", "3"]

# hand cases: lists of (a, b); None for a missing cell
FIXED = [
    [(1.5, 1.0), (2.1, 2.5), (3.7, 2.0), (0.3, 1.1), (5.5, 3.0), (2.2, 2.9), (4.1, 3.0), (6.0, 1.0)],
    [(1.0, 1.0), (2.0, 2.0), (3.0, 2.0), (4.0, 3.0), (5.0, 5.0), (5.0, 6.0), (6.0, 6.0), (7.0, 9.0), (8.0, 2.0)],
    [(1.0, 2.0)],
    [(2.0, 1.0)],
    [(3.0, 3.0), (4.0, 4.0)],
    [(1.0, 2.0), (3.0, 1.0)],
    [(1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)],
    [(5e-324, 0.0), (1e-323, 5e-324), (0.0, 0.0), (2e-323, 5e-324)],
    [(1e308, -1e308), (1.0, 2.0)],
    [(-1e308, 1e308), (1.0, 2.0)],
    [(1e300, 1e299), (-1e300, 3e299), (2e300, 1.0)],
    [(1.0, 1.0000000000000002), (2.0, 2.0000000000000004), (3.0, 2.9999999999999996), (4.0, 4.0)],
    [(float(i), float(i % 7)) for i in range(1, 52)],
    [(float(i) + 0.5, float((i * 3) % 11)) for i in range(1, 51)],
    [(float(i) + 0.5, float((i * 3) % 11)) for i in range(1, 52)],
    [(0.1 * i, 0.3) for i in range(1, 12)],
    [(1.0, None), (None, 2.0), (3.0, 4.0), (4.0, 7.0), (5.0, 5.5), (2.0, 1.0)],
]


def run(binary, tmp, names, rows, extra):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "test", "wilcoxon", "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=600)


class Stat:
    def __init__(self):
        self.ctl = {}
        self.p = {"max": 0.0, "n": 0}
        self.exact = {"n": 0, "differs": 0}
        self.sp = {"w": 0.0, "pe": 0.0, "pn": 0.0, "n": 0}
        self.checked = 0
        self.refused = 0
        self.reasons = {}

    def control(self, name, outside):
        c = self.ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0


def gen_cols(r, n, kind):
    """one column of n cells of the kind: 'lattice' (small integers, many equal and zero differences), 'uniform', 'ints', 'offset', 'wide', 'tiny', 'gappy', 'const'"""
    out = []
    for i in range(n):
        if r.next() % 30 == 0:
            out.append("")
            continue
        if kind == "lattice":
            out.append(str(int(r.next() % 6)))
        elif kind == "ints":
            out.append(str(int(r.next() % 41) - 20))
        elif kind == "offset":
            out.append(repr(1e8 + r.u()))
        elif kind == "wide":
            out.append(repr(math.ldexp(1.0 + r.u(), int(r.next() % 600) - 300) * (1 if r.next() & 1 else -1)))
        elif kind == "tiny":
            out.append(repr(math.ldexp(1.0 + r.u(), -1000 - int(r.next() % 60))))
        elif kind == "const":
            out.append("0.1")
        else:
            out.append(repr(round(r.u() * 10, 3)))
    return out


def gen(r, fi):
    if fi < len(FIXED):
        rows = [[fmt(a), fmt(b)] for a, b in FIXED[fi]]
        return ["a", "b"], rows, ["fixed", "fixed"]
    n = [1, 2, 3, 5, 8, 12, 20, 40, 50, 51, 60, 200][int(r.next() % 12)]
    kinds = ["lattice", "lattice", "uniform", "uniform", "ints", "offset", "wide", "tiny", "const"]
    ka = kinds[int(r.next() % len(kinds))]
    shared = r.next() % 2 == 0
    ca = gen_cols(r, n, ka)
    if shared:
        # b = a plus a small change (often none): zero differences, ties, and a shared part
        cb = []
        for v in ca:
            if v == "" or r.next() % 25 == 0:
                cb.append("")
            elif r.next() % 4 == 0:
                cb.append(v)
            else:
                try:
                    cb.append(repr(float(v) + (int(r.next() % 9) - 3) * 0.5))
                except ValueError:
                    cb.append(v)
        kb = "shared"
    else:
        kb = kinds[int(r.next() % len(kinds))]
        cb = gen_cols(r, n, kb)
    names = ["a_" + ka, "b_" + kb]
    rows = [[ca[i], cb[i]] for i in range(n)]
    return names, rows, [ka, kb if not shared else "shared"]


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
    for n in range(1, 9):
        if ot.wx_counts(n) != ot.wx_brute(n):
            print("gate_wilcoxon: the oracle's counts differ from the enumeration at n = %d" % n)
            return 2
    files = 0
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        nf = len(FIXED)
        for fi in range(a.files):
            names, rows, shapes = gen(r, fi % nf if fi < 9 * nf else fi)
            n = len(rows)
            if not column_numeric(rows, 0) or not column_numeric(rows, 1):
                continue
            fixed = fi < 9 * nf
            onesample = r.next() % 3 == 0 and not fixed
            # the hand cases run under all three methods for zeros and all three exact modes
            zeros = ["wilcox", "pratt", "zsplit"][(fi // nf) % 3] if fixed else r.pick(["wilcox", "wilcox", "pratt", "zsplit"])
            exact = ["auto", "yes", "no"][(fi // (3 * nf)) % 3] if fixed else r.pick(["auto", "auto", "no", "yes"])
            cont = r.pick(["auto", "auto", "none"])
            mu0 = r.pick(MU0S)
            if onesample:
                extra = ["--cols", names[0], "--mu0", mu0]
                pairs_in = [(parse(rw[0]), float(mu0)) for rw in rows if parse(rw[0]) is not None]
            else:
                extra = ["--pairs", "%s:%s" % (names[0], names[1])]
                pairs_in = [(parse(rw[0]), parse(rw[1])) for rw in rows if parse(rw[0]) is not None and parse(rw[1]) is not None]
            extra += ["--zeros", zeros, "--exact", exact, "--continuity", cont]
            correction = r.pick(["none", "none", "holm"])
            if correction != "none":
                extra += ["--correct", correction]
            w = ot.wilcoxon(pairs_in, zeros, exact, cont == "auto")
            try:
                p = run(a.bin, tmp, names, rows, extra)
            except subprocess.TimeoutExpired:
                print("gate_wilcoxon: the program did not finish (hang)")
                return 2
            if exact == "yes" and w["reason"] is None and not w["exact_available"]:
                st.refused += 1
                doc = json.loads(p.stdout) if p.stdout else {}
                if p.returncode != 2 or doc.get("error", {}).get("rule") != "stats.exact-unavailable":
                    bad.append((fi, "--exact yes where the exact p-value is not available: want exit 2 stats.exact-unavailable, got exit %d %r" % (p.returncode, p.stdout[:160])))
                continue
            if p.returncode != 0:
                keep = pathlib.Path(tempfile.gettempdir()) / ("gate_wilcoxon_failed_%d.csv" % fi)
                keep.write_text((tmp / "f.csv").read_text())
                print("gate_wilcoxon: exit %d on file %d: %s (the file is %s, args %r)" % (p.returncode, fi, p.stdout[:240], keep, extra))
                return 2
            files += 1
            doc = json.loads(p.stdout)["data"]
            g = doc["tests"][0]
            tag = "wilcoxon %s %s %s %s" % ("1s" if onesample else "pairs", zeros, exact, cont)
            st.checked += 1
            st.reasons[w["reason"]] = st.reasons.get(w["reason"], 0) + 1
            if doc.get("test") != "wilcoxon" or doc.get("rows") != n or doc.get("adjustment") != correction:
                bad.append((fi, "the header: %r" % {k: doc.get(k) for k in ("test", "rows", "adjustment")}))
            if correction != "none":
                want_adj = ot.adjust([g["p_value"]], correction)
                fam = doc.get("family") or {}
                if fam.get("correction") != correction or not all(gd.same(a_, b_) for a_, b_ in zip(fam.get("p_adjusted", []), want_adj)):
                    bad.append((fi, "the family: want %s %r got %r" % (correction, want_adj, fam)))
            if g.get("reason") != w["reason"]:
                bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
            if (g.get("column") != names[0]) if onesample else (g.get("pair") != [names[0], names[1]]):
                bad.append((fi, "%s: the test names %r" % (tag, g.get("column", g.get("pair")))))
            if "level" in doc or (doc.get("mu0") is None) != (not onesample) or (onesample and not gd.same(doc.get("mu0"), float(mu0))):
                bad.append((fi, "%s: the header has level %r mu0 %r" % (tag, doc.get("level"), doc.get("mu0"))))
            if g["zero_method"] != zeros or g.get("continuity") != (cont == "auto") or g.get("method") != "wilcoxon":
                bad.append((fi, "%s echo: %r" % (tag, {k: g.get(k) for k in ("zero_method", "continuity", "method")})))
            if w["reason"] in ("too-few-rows", "difference-beyond-double"):
                if g.get("statistic") is not None or g.get("p_value") is not None:
                    bad.append((fi, "%s: want nulls got %r" % (tag, g)))
                if w["reason"] == "too-few-rows" and g["n"]["used"] != 0:
                    bad.append((fi, "%s used: %r" % (tag, g["n"])))
                continue
            if g["n"] != {"used": w["n"], "nonzero": w["nz"], "zeros": w["zeros"]}:
                bad.append((fi, "%s n: want %r got %r" % (tag, {"used": w["n"], "nonzero": w["nz"], "zeros": w["zeros"]}, g["n"])))
            if True:
                if not gd.same(g["statistic"], w["wplus"]) or not gd.same(g["w"]["plus"], w["wplus"]) or not gd.same(g["w"]["minus"], w["wminus"]):
                    bad.append((fi, "%s W: want %r %r got %r %r" % (tag, w["wplus"], w["wminus"], g.get("statistic"), g.get("w"))))
                if g["ties"] != w["ties"]:
                    bad.append((fi, "%s ties: want %r got %r" % (tag, w["ties"], g["ties"])))
                es = g.get("effect_size")
                if w["rb"] is None:
                    if es is not None:
                        bad.append((fi, "%s effect_size: want null got %r" % (tag, es)))
                elif es is None or not gd.same(es.get("value"), w["rb"]):
                    bad.append((fi, "%s rank_biserial: want %r got %r" % (tag, w["rb"], es)))
            if w["z"] is None:
                if g.get("z") is not None:
                    bad.append((fi, "%s z: want null got %r" % (tag, g.get("z"))))
            elif not gd.same(g.get("z"), w["z"]):
                bad.append((fi, "%s z: want %r got %r" % (tag, w["z"], g.get("z"))))
            if w["p"] is None:
                if g.get("p_value") is not None:
                    bad.append((fi, "%s p_value: want null got %r" % (tag, g.get("p_value"))))
            else:
                if g.get("p_method") != w["p_method"]:
                    bad.append((fi, "%s p_method: want %r got %r" % (tag, w["p_method"], g.get("p_method"))))
                if g.get("p_value") is None:
                    bad.append((fi, "%s p_value: want %s got null" % (tag, mpmath.nstr(w["p"], 17))))
                elif w["p_method"] == "exact":
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
                if w["z2"] is not None and w["p"] > mpmath.mpf("1e-250") and w["wplus"] is not None:
                    if w["p_method"] == "normal":
                        alt = ot.wilcoxon(pairs_in, zeros, "no", cont != "auto")
                        if alt["p"] is not None:
                            st.control("other-continuity", not ot.p_ok(float(alt["p"]), w["p"], 1.0))
                        if w["zeros"] > 0 and zeros != "wilcox":
                            alt = ot.wilcoxon(pairs_in, "wilcox", "no", cont == "auto")
                            if alt["p"] is not None:
                                st.control("other-zero-method", not ot.p_ok(float(alt["p"]), w["p"], 1.0))
                    if w["p_method"] == "exact":
                        alt = ot.wilcoxon(pairs_in, zeros, "no", cont == "auto")
                        if alt["p"] is not None:
                            st.control("normal-for-exact", not ot.p_ok(float(alt["p"]), w["p"], 1.0))
                    st.control("one-sided", not ot.p_ok(float(w["p"] / 2), w["p"], 1.0))
                if a.scipy and w["reason"] is None and w["wplus"] is not None and w["nz"] > 0 and all(Fraction(x) - Fraction(y) == Fraction(x - y) for x, y in pairs_in):
                    import scipy.stats as ss
                    xs = [x for x, _ in pairs_in]
                    ys = [y for _, y in pairs_in]
                    try:
                        res = ss.wilcoxon(xs, ys, zero_method=zeros, method="exact" if w["p_method"] == "exact" else "approx", correction=(cont == "auto"))
                    except Exception:
                        res = None
                    if res is not None:
                        want_stat = min(float(w["wplus"]), float(w["wminus"])) if zeros != "zsplit" or True else None
                        du = abs(res.statistic - want_stat) / max(abs(want_stat), 1.0)
                        dp = abs(res.pvalue - g["p_value"]) / max(abs(res.pvalue), 1e-300)
                        st.sp["n"] += 1
                        st.sp["w"] = max(st.sp["w"], du)
                        key = "pe" if w["p_method"] == "exact" else "pn"
                        st.sp[key] = max(st.sp[key], dp)
                        if du > 1e-12 or dp > (1e-12 if key == "pe" else 1e-9):
                            bad.append((fi, "%s scipy: dW %g dp %g (scipy %r %r, ours %r %r)" % (tag, du, dp, res.statistic, res.pvalue, want_stat, g["p_value"])))
            if any(b[0] == fi for b in bad):
                pathlib.Path(tempfile.gettempdir(), "gate_wilcoxon_bad_%d.csv" % fi).write_text(",".join(names) + "\n" + "\n".join(",".join(rw) for rw in rows) + "\n")
                pathlib.Path(tempfile.gettempdir(), "gate_wilcoxon_bad_%d.args" % fi).write_text(" ".join(extra))
            if n > 1 and w["reason"] != "difference-beyond-double":
                p2 = run(a.bin, tmp, names, shuffled(r, rows), extra)
                if p2.stdout != p.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
                if not onesample:
                    # the pair named the other way round and the values times 8
                    p3 = run(a.bin, tmp, names, rows, ["--pairs", "%s:%s" % (names[1], names[0])] + extra[2:])
                    if p3.returncode == 0:
                        g3 = json.loads(p3.stdout)["data"]["tests"][0]
                        neg = lambda x: None if x is None else (-x if x != 0 else 0.0)
                        ok = g3.get("p_value") == g.get("p_value") and g3.get("ties") == g.get("ties") and g3.get("n") == g.get("n")
                        if g.get("w") is not None and zeros != "zsplit":
                            ok = ok and g3["w"]["plus"] == g["w"]["minus"] and g3["w"]["minus"] == g["w"]["plus"]
                        ok = ok and gd.same(neg(g3.get("z")), neg(neg(g.get("z")))) is not None
                        if g.get("z") is not None:
                            ok = ok and gd.same(g3["z"], neg(g["z"]))
                        if g.get("effect_size") is not None:
                            ok = ok and gd.same(g3["effect_size"]["value"], neg(g["effect_size"]["value"]))
                        if not ok:
                            bad.append((fi, "the pair named the other way round does not mirror the answer: %r against %r" % (g, g3)))
                    trows = [[(repr(parse(c) * 8.0) if c != "" else c) for c in rw] for rw in rows]
                    p5 = run(a.bin, tmp, names, trows, extra if not onesample else ["--cols", names[0], "--mu0", repr(float(mu0) * 8.0)] + extra[4:])
                    if p5.returncode == 0:
                        g5 = json.loads(p5.stdout)["data"]["tests"][0]
                        keys = ("statistic", "w", "z", "p_value", "p_method", "effect_size", "n", "ties", "reason")
                        if any(g.get(k_) != g5.get(k_) for k_ in keys) and not any(s_ in ("wide", "tiny") for s_ in shapes):
                            bad.append((fi, "the values times 8 change the answer: %r against %r" % (g, g5)))
    for k, (outside, total) in sorted(st.ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_wilcoxon: %d files, %d tests compared (%s), %d refusals expected and seen, %d differ; exact p-values equal to the oracle's rational rounded once: %d of %d; normal p-values against mpmath: largest relative error %.3g over %d; controls outside the bound: %s%s" % (
        files, st.checked, ", ".join("%s %d" % (k or "answered", v) for k, v in sorted(st.reasons.items(), key=lambda x: str(x[0]))), st.refused, len(bad), st.exact["n"] - st.exact["differs"], st.exact["n"], st.p["max"], st.p["n"],
        ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(st.ctl.items())),
        ("; scipy on %d tests: max relative dW %.3g, exact dp %.3g, normal dp %.3g" % (st.sp["n"], st.sp["w"], st.sp["pe"], st.sp["pn"])) if a.scipy else ""))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
