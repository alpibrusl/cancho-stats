#!/usr/bin/env python3
"""S5.1 gates I-A and the invariances (docs/resample.md section 6) for `stats ci prop`:

    python3 scripts/gate_ci_prop.py [--bin build/stats] [--files 200] [--seed 1]

Generated CSV files (an outcome column with a success label among others, a group column of two labels, empty cells in both, trials from 0 to 3,000, proportions of 0, 1 and everything between, levels from 0.5 to 0.999 and one of nine
decimals) run with `ci prop`, one sample (`--method wilson`, the default, and `clopper-pearson`) and two samples (`--group`, Newcombe). The Wilson and Newcombe `interval` objects must be BYTE EQUAL to the `interval` of `stats test prop`
on the same file (the gate of the test is `gate_test4.py`, against mpmath), the counts equal the oracle's. Clopper-Pearson against scripts/oracle_test.py's `cp_interval`: the quantiles of Beta(x, n - x + 1) and the upper tail of
Beta(x + 1, n - x) at the DOUBLE alpha/2 the program uses, by bisection on the regularised incomplete beta function of mpmath at 50 digits; each endpoint within the bound of docs/numerics.md 9.2, relative
1e-13 + 1e-15 |ln(alpha/2)| (kappa = 1: both shapes are at least 1); x = 0 has the lower endpoint 0 and x = n the upper endpoint 1, exactly. Also: the rows shuffled give the same bytes; two columns of one file give the entries of the two
files; `--method newcombe` without `--group` and `--method clopper-pearson` with it are refused (args.bad-value). CONTROLS from the same data must leave the bound on at least a tenth of the cases: the Wald interval, Wilson's in place of
Clopper-Pearson's, and the normal-quantile lower endpoint.

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
from gate_anova import shuffled  # noqa: E402

LEVELS = ["0.95", "0.95", "0.9", "0.99", "0.5", "0.999", "0.8", "0.123456789"]
OUTCOMES = [("yes", ["no"]), ("1", ["0"]), ("a b", ["c", "d"]), ("x,y", ["z"]), ("T", ["F", "U"])]


def run(binary, tmp, names, rows, extra, verb=("ci", "prop")):
    text = ",".join(names) + "\n" + "\n".join(",".join('"%s"' % c if "," in c else c for c in r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, verb[0], verb[1], "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=300)


def gen(r, fi):
    n = [0, 1, 2, 3, 5, 10, 20, 50, 200, 1000, 3000][int(r.next() % 11)] if fi >= 12 else [1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20][fi]
    succ, others = OUTCOMES[int(r.next() % len(OUTCOMES))]
    p = [0.0, 1.0, 0.5, 0.1, 0.9, 0.02, 0.98, 0.3][int(r.next() % 8)]
    labels = sorted(["g1", "g2"])
    rows = []
    for _ in range(n):
        y = "" if r.u() < 0.04 else (succ if r.u() < p else others[int(r.next() % len(others))])
        g = "" if r.u() < 0.04 else labels[int(r.next() % 2)]
        rows.append([y, g, "" if r.u() < 0.04 else (succ if r.u() < 0.5 else others[0])])
    if n >= 2:
        rows[0][1] = labels[0]
        rows[1][1] = labels[1]
    return ["y", "g", "w"], rows, succ, labels


def cp_interval(x, n, a2):
    """the Clopper-Pearson endpoints at the double a2 = alpha / 2: (low, high) as mpf; x = 0 gives low 0 and x = n high 1"""
    mp = mpmath.mp
    mp.dps = 50
    a2 = mpmath.mpf(a2)

    def bisect(f, lo, hi):
        # f increasing in x from f(lo) < 0 to f(hi) > 0
        for _ in range(220):
            mid = (lo + hi) / 2
            if f(mid) < 0:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    low = mpmath.mpf(0)
    high = mpmath.mpf(1)
    if x > 0:
        # I_t(x, n - x + 1) = a2
        low = bisect(lambda t: mpmath.betainc(x, n - x + 1, 0, t, regularized=True) - a2, mpmath.mpf(0), mpmath.mpf(1))
    if x < n:
        # the upper tail 1 - I_t(x + 1, n - x) = a2, i.e. I_t(x + 1, n - x) = 1 - a2
        high = bisect(lambda t: mpmath.betainc(x + 1, n - x, 0, t, regularized=True) - (1 - a2), mpmath.mpf(0), mpmath.mpf(1))
    return low, high


def half_alpha(level):
    num, k = ot.level_parts(level)
    scale = 10 ** k
    return (scale - num) / (2 * scale)   # the double, as the program computes float(scale - num) / float(2 scale)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    ctl = {}
    files = tests = cp_n = 0
    worst = 0.0

    def control(name, outside):
        c = ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0

    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            names, rows, succ, labels = gen(r, fi)
            level = LEVELS[int(r.next() % len(LEVELS))]
            n = len(rows)
            ys = [rw[0] for rw in rows if rw[0] != ""]
            x = sum(1 for v in ys if v == succ)
            nn = len(ys)
            base = ["--cols", "y,w", "--success", succ, "--level", level]
            # Wilson: the interval of `test prop` byte for byte
            pw = run(a.bin, tmp, names, rows, base)
            pt = run(a.bin, tmp, names, rows, base, verb=("test", "prop"))
            if pw.returncode != 0 or pt.returncode != 0:
                print("gate_ci_prop: exit %d / %d on file %d: %s %s" % (pw.returncode, pt.returncode, fi, pw.stdout[:200], pt.stdout[:200]))
                return 2
            files += 1
            gw = json.loads(pw.stdout)["data"]["intervals"]
            gt = json.loads(pt.stdout)["data"]["tests"]
            for k, (g, h) in enumerate(zip(gw, gt)):
                tests += 1
                if g["interval"] != h["interval"]:
                    bad.append((fi, "wilson %s: ci %r against test %r" % (names[[0, 2][k]], g["interval"], h["interval"])))
                if g["n"] != h["n"]:
                    bad.append((fi, "wilson n: %r against %r" % (g["n"], h["n"])))
                if g.get("reason") != h.get("reason") and not (g.get("reason") is None and h.get("reason") == "variance-zero"):
                    bad.append((fi, "wilson reason: %r against %r" % (g.get("reason"), h.get("reason"))))
            # Clopper-Pearson, first column y
            pc = run(a.bin, tmp, names, rows, base + ["--method", "clopper-pearson"])
            if pc.returncode != 0:
                bad.append((fi, "clopper-pearson: exit %d %r" % (pc.returncode, pc.stdout[:160])))
                continue
            gc = json.loads(pc.stdout)["data"]["intervals"][0]
            if gc["n"]["n"] != nn or gc["n"]["successes"] != x or gc["n"]["missing"] != n - nn:
                bad.append((fi, "cp counts: %r want n %d x %d missing %d" % (gc["n"], nn, x, n - nn)))
            if nn == 0:
                if gc.get("interval") is not None or gc.get("reason") != "too-few-rows":
                    bad.append((fi, "cp with no trial: %r" % gc))
            else:
                cp_n += 1
                a2 = half_alpha(level)
                lo, hi = cp_interval(x, nn, a2)
                iv = gc["interval"]
                tol_rel = mpmath.mpf("1e-13") + mpmath.mpf("1e-15") * abs(mpmath.log(mpmath.mpf(a2)))
                for key, want in (("low", lo), ("high", hi)):
                    got = mpmath.mpf(iv[key])
                    if want == 0 or want == 1:
                        ok = got == want
                    else:
                        ok = abs(got - want) <= tol_rel * abs(want)
                        worst = max(worst, float(abs(got - want) / abs(want) / tol_rel))
                    if not ok:
                        bad.append((fi, "cp %s x %d n %d level %s: want %s got %r" % (key, x, nn, level, mpmath.nstr(want, 17), iv[key])))
                if iv["method"] != "clopper-pearson" or iv["of"] != "proportion" or not gd.same(iv["estimate"], x / nn):
                    bad.append((fi, "cp object: %r" % iv))
                # controls
                q = abs(float(mpmath.sqrt(2) * mpmath.erfinv(1 - 2 * a2)))
                ph = x / nn
                wald_lo, wald_hi = ph - q * math.sqrt(ph * (1 - ph) / nn), ph + q * math.sqrt(ph * (1 - ph) / nn)
                wl = json.loads(pw.stdout)["data"]["intervals"][0]["interval"]
                if 0 < x < nn:
                    def far(got, want):
                        return abs(got - float(want)) > float(tol_rel) * abs(float(want))
                    control("wald", far(wald_lo, lo) or far(wald_hi, hi))
                    control("wilson-for-cp", far(wl["low"], lo) or far(wl["high"], hi))
            # the rows shuffled, the columns in another order
            if n > 1:
                p2 = run(a.bin, tmp, names, shuffled(r, rows), base + ["--method", "clopper-pearson"])
                if p2.stdout != pc.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
                p3 = run(a.bin, tmp, names, rows, ["--cols", "w,y", "--success", succ, "--level", level, "--method", "clopper-pearson"])
                d3 = json.loads(p3.stdout)["data"]["intervals"]
                d0 = json.loads(pc.stdout)["data"]["intervals"]
                if [d3[0], d3[1]] != [dict(d0[1]), dict(d0[0])] and not (d3[0]["column"] == "w" and d3[0]["interval"] == d0[1]["interval"] and d3[1]["interval"] == d0[0]["interval"]):
                    bad.append((fi, "--cols in another order changes an entry"))
            # two samples
            if n >= 2:
                two = ["--cols", "y,w", "--success", succ, "--group", "g", "--level", level]
                p4 = run(a.bin, tmp, names, rows, two)
                p5 = run(a.bin, tmp, names, rows, two, verb=("test", "prop"))
                glev = sorted(set(rw[1] for rw in rows if rw[1] != ""))
                if len(glev) == 2:
                    if p4.returncode != 0 or p5.returncode != 0:
                        bad.append((fi, "newcombe: exit %d / %d" % (p4.returncode, p5.returncode)))
                    else:
                        for k, (g, h) in enumerate(zip(json.loads(p4.stdout)["data"]["intervals"], json.loads(p5.stdout)["data"]["tests"])):
                            tests += 1
                            if g["interval"] != h["interval"] or g["n"] != h["n"] or g["levels"] != h["levels"]:
                                bad.append((fi, "newcombe %d: ci %r against test %r" % (k, g, h)))
        # refusals
        names, rows, succ, labels = gen(SM64(9), 14)
        for extra, what in ((["--cols", "y", "--success", succ, "--method", "newcombe"], "newcombe without --group"),
                            (["--cols", "y", "--success", succ, "--group", "g", "--method", "clopper-pearson"], "clopper-pearson with --group"),
                            (["--cols", "y", "--success", succ, "--group", "g", "--method", "wilson"], "wilson with --group"),
                            (["--cols", "y", "--success", succ, "--method", "pearson"], "pearson for a proportion")):
            p = run(a.bin, tmp, names, rows, extra)
            doc = json.loads(p.stdout) if p.stdout else {}
            if p.returncode != 2 or doc.get("error", {}).get("rule") != "args.bad-value":
                bad.append((-2, "%s: want exit 2 args.bad-value, got exit %d %r" % (what, p.returncode, p.stdout[:120])))
    for k, (outside, total) in sorted(ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_ci_prop: %d files, %d entries compared with `test prop`, %d Clopper-Pearson intervals against mpmath (largest error %.3g of the bound), %d problems; controls outside the bound: %s" % (
        files, tests, cp_n, worst, len(bad), ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(ctl.items()))))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
