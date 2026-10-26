#!/usr/bin/env python3
"""S4.7 gates for `stats corr --method kendall --exact auto|yes|no` (docs/rank.md section 15):

    python3 scripts/gate_corr_exact.py [--bin build/stats] [--files 200] [--seed 1] [--scipy]

Generated CSV files (two to four columns: permutations of distinct values with a planted dependence, columns with ties, missing cells, sizes from 2 to 100 so that both sides of the limit of 60 and of the bound n >= 3 are met)
against scripts/oracle_test.py, which takes S from the pairs (no sorting), the tie groups by counting, counts the permutations of n by inversions with the product of the polynomials 1 + q + .. + q^(i - 1) (checked first against
the enumeration of every permutation for n <= 8), and rounds the exact two-sided p-value once. For each pair and each of `--exact` absent, `no`, `auto`, `yes`: with `no` or absent the bytes of the answer are the same as without the
flag (nothing of S2 changes); with `auto` a pair with no ties in either column and 3 to 60 rows has `p_method: "exact"` and the p-value EQUAL to the oracle's rational rounded once, every other pair the bytes of the default run;
with `yes` the same where every pair can, and otherwise (any tie, or more than 60 rows, in a pair) a refusal, stats.exact-unavailable, exit 2; `--exact yes` with pearson or spearman is that refusal too. CONTROLS from the same
data must leave the bound on at least a tenth of the cases: the normal-approximation p-value, and half the exact p-value. With --scipy: `kendalltau(method="exact")` agrees to 1e-12 relative where the data have no ties.

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

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
import gate_describe as gd  # noqa: E402
import oracle_test as ot  # noqa: E402
from gate_test import parse  # noqa: E402

FIXED_N = [2, 3, 3, 4, 5, 6, 7, 10, 25, 59, 60, 61, 100]


def run(binary, tmp, names, rows, extra=(), method="kendall"):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "corr", "--method", method, "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=300)


def gen(r, fi):
    n = FIXED_N[fi] if fi < len(FIXED_N) else [2, 3, 4, 5, 8, 12, 20, 40, 60, 61, 80][int(r.next() % 11)]
    ncol = 2 + int(r.next() % 3)
    cols = []
    base = list(range(n))
    for j in range(ncol):
        kind = "ties" if (r.next() % 4 == 0 and fi >= len(FIXED_N)) else "perm"
        # a permutation correlated with the base order: sort by base + noise
        if kind == "perm":
            noise = [(i + (r.u() - 0.5) * n * [0.0, 0.5, 2.0, 20.0][int(r.next() % 4)], i) for i in base]
            order = sorted(noise)
            vals = [0] * n
            for rank, (_, i) in enumerate(order):
                vals[i] = rank
            cols.append([str(v * 3 - 7) for v in vals])
        else:
            cols.append([str(int(r.next() % 4)) for _ in base])
    # missing cells
    for c in cols:
        for i in range(n):
            if r.u() < 0.03 and fi >= len(FIXED_N):
                c[i] = ""
    names = ["c%d" % j for j in range(ncol)]
    rows = [[cols[j][i] for j in range(ncol)] for i in range(n)]
    return names, rows


def pair_truth(rows, ja, jb):
    pts = [(parse(rw[ja]), parse(rw[jb])) for rw in rows if rw[ja] != "" and rw[jb] != ""]
    n = len(pts)
    s = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = (Fraction(pts[i][0]) > Fraction(pts[j][0])) - (Fraction(pts[i][0]) < Fraction(pts[j][0]))
            dy = (Fraction(pts[i][1]) > Fraction(pts[j][1])) - (Fraction(pts[i][1]) < Fraction(pts[j][1]))
            s += dx * dy
    def groups(vals):
        v = [Fraction(x) for x in vals]
        return sum(1 for x in set(v) if v.count(x) > 1)
    return n, s, groups([p[0] for p in pts]), groups([p[1] for p in pts]), pts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scipy", action="store_true")
    a = ap.parse_args()
    r = SM64(a.seed)
    for n in range(1, 9):
        if ot.mahonian(n) != ot.kendall_brute(n):
            print("gate_corr_exact: the oracle's counts differ from the enumeration at n = %d" % n)
            return 2
    bad = []
    ctl = {}
    files = pairs = exact_pairs = exact_differs = refusals = 0
    sp = {"p": 0.0, "n": 0}

    def control(name, outside):
        c = ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0

    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            names, rows = gen(r, fi)
            ncol = len(names)
            base = run(a.bin, tmp, names, rows)
            if base.returncode != 0:
                print("gate_corr_exact: the default run failed on file %d: %s" % (fi, base.stdout[:200]))
                return 2
            files += 1
            bp = json.loads(base.stdout)["data"]["pairs"]
            # the flag absent and "no" give the same bytes
            pno = run(a.bin, tmp, names, rows, ["--exact", "no"])
            if pno.stdout != base.stdout:
                bad.append((fi, "--exact no changes the bytes of the default run"))
            truths = {}
            k = 0
            can = True
            for ja in range(ncol):
                for jb in range(ja + 1, ncol):
                    n, s, gx, gy, pts = pair_truth(rows, ja, jb)
                    truths[k] = (n, s, gx, gy, pts)
                    if not (gx == 0 and gy == 0 and n <= 60):
                        can = False
                    k += 1
            pauto = run(a.bin, tmp, names, rows, ["--exact", "auto"])
            if pauto.returncode != 0:
                bad.append((fi, "--exact auto: exit %d %r" % (pauto.returncode, pauto.stdout[:160])))
                continue
            ap_ = json.loads(pauto.stdout)["data"]["pairs"]
            for k, (g, g0) in enumerate(zip(ap_, bp)):
                n, s, gx, gy, pts = truths[k]
                pairs += 1
                if gx == 0 and gy == 0 and 3 <= n <= 60 and g0.get("p_value") is not None:
                    exact_pairs += 1
                    want = ot.rounded(ot.kendall_exact(n, s))
                    if g.get("p_method") != "exact" or not gd.same(g.get("p_value"), want):
                        exact_differs += 1
                        bad.append((fi, "pair %d auto: want exact %r got %r %r" % (k, want, g.get("p_method"), g.get("p_value"))))
                    # everything else is the default run's
                    for key in ("x", "y", "n", "tau", "statistic", "s", "ties"):
                        if g.get(key) != g0.get(key):
                            bad.append((fi, "pair %d auto: %s differs from the default run: %r against %r" % (k, key, g.get(key), g0.get(key))))
                    control("normal-approximation", not gd.same(g0["p_value"], want) and not (abs(g0["p_value"] - want) <= 1e-12 * want))
                    control("half-p", not gd.same(want / 2, want))
                    if a.scipy:
                        import scipy.stats as ss
                        res = ss.kendalltau([p[0] for p in pts], [p[1] for p in pts], method="exact", variant="b")
                        dp = abs(res.pvalue - g["p_value"]) / max(abs(res.pvalue), 1e-300)
                        sp["n"] += 1
                        sp["p"] = max(sp["p"], dp)
                        if dp > 1e-12:
                            bad.append((fi, "pair %d scipy: dp %g" % (k, dp)))
                elif g != g0:
                    bad.append((fi, "pair %d auto: not exact-eligible, want the default bytes: %r against %r" % (k, g, g0)))
            pyes = run(a.bin, tmp, names, rows, ["--exact", "yes"])
            if can:
                if pyes.stdout != pauto.stdout:
                    bad.append((fi, "--exact yes where every pair can: want the bytes of auto"))
            else:
                refusals += 1
                doc = json.loads(pyes.stdout) if pyes.stdout else {}
                if pyes.returncode != 2 or doc.get("error", {}).get("rule") != "stats.exact-unavailable":
                    bad.append((fi, "--exact yes where a pair cannot: want exit 2 stats.exact-unavailable, got exit %d %r" % (pyes.returncode, pyes.stdout[:160])))
        # pearson and spearman with yes
        names, rows = gen(SM64(5), 5)
        for m in ("pearson", "spearman"):
            p = run(a.bin, tmp, names, rows, ["--exact", "yes"], method=m)
            doc = json.loads(p.stdout) if p.stdout else {}
            refusals += 1
            if p.returncode != 2 or doc.get("error", {}).get("rule") != "stats.exact-unavailable":
                bad.append((-2, "--exact yes with %s: want exit 2 stats.exact-unavailable, got exit %d" % (m, p.returncode)))
            for mode in ("auto", "no"):
                q = run(a.bin, tmp, names, rows, ["--exact", mode], method=m)
                q0 = run(a.bin, tmp, names, rows, [], method=m)
                if q.stdout != q0.stdout:
                    bad.append((-2, "--exact %s changes %s" % (mode, m)))
    for k, (outside, total) in sorted(ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_corr_exact: %d files, %d pairs compared, %d exact pairs (%d differ), %d refusals expected and seen, %d problems; controls outside the bound: %s%s" % (
        files, pairs, exact_pairs, exact_differs, refusals, len(bad), ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(ctl.items())),
        ("; scipy on %d pairs: max relative dp %.3g" % (sp["n"], sp["p"])) if a.scipy else ""))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
