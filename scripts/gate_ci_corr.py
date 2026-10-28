#!/usr/bin/env python3
"""S5.1c gate I-A (docs/resample.md sections 3 and 11) for `stats ci corr`:

    python3 scripts/gate_ci_corr.py [--bin build/stats] [--files 200] [--seed 1]

The files of scripts/gate_corr.py (the numeric shapes of gate_describe.py, correlated pairs, ties, missing cells, huge and tiny values, constants, sizes from 1 to 200) run with `ci corr` in both `--missing` modes and levels from 0.5 to
0.999. For every pair of numeric columns: n, the estimate (the double r of `corr`, bit for bit from the oracle's exact N, Lx, Ly), the reason (`constant-column`, `perfect-correlation` where |r| = 1 is decided by the integers, `too-few-rows`
for n < 4) and the Fisher interval against mpmath at 60 digits FROM THE DOUBLE r the program prints: z = atanh r, standard error 1 / sqrt(n - 3), q the normal quantile of (1 + level) / 2 at 60 digits, endpoints tanh(z -+ q se). The
bound is fixed before the first run: relative 1e-12 of max(|endpoint|, (1 - endpoint^2) max(|z|, q se)) -- the second term is the rounding of the argument z -+ q se, which a tanh near 0 or near 1 turns into a larger relative error of
the endpoint, and the normal quantile's own 2e-13 enters through q se. The interval depends on the DOUBLE r (a conditioning of 1 / (1 - r^2) near |r| = 1, stated in docs/resample.md section 11). Also: the rows shuffled give the same
bytes; `--method spearman` and `kendall` are refused (args.bad-value). CONTROLS from the same data must leave the bound on at least a tenth of the cases: n - 2 in place of n - 3, and the endpoints of z without tanh.

Exit 0: every number is the expected one. Exit 1: a difference. Exit 2: the harness failed.
"""
import argparse
import json
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
import gate_corr as gc  # noqa: E402
import oracle_test as ot  # noqa: E402
from gate_anova import shuffled  # noqa: E402

LEVELS = ["0.95", "0.95", "0.9", "0.99", "0.5", "0.999", "0.8", "0.123456789"]
mpmath.mp.dps = 60


def run(binary, tmp, names, rows, extra):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "ci", "corr", "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=300)


def fisher(r, n, q, df_off=3):
    z = mpmath.atanh(mpmath.mpf(r))
    se = 1 / mpmath.sqrt(n - df_off)
    return z, q * se, mpmath.tanh(z - q * se), mpmath.tanh(z + q * se)


def within(got, want, z, qse):
    scale = max(abs(want), (1 - want * want) * max(abs(z), qse))
    return abs(mpmath.mpf(got) - want) <= mpmath.mpf("1e-12") * scale


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    ctl = {}
    pairs = intervals = 0
    worst = 0.0

    def control(name, outside):
        c = ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0

    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            if fi < len(gc.FIXED):
                cols = gc.FIXED[fi]
                shapes = ["fixed"] * len(cols)
                ncols = len(cols)
                n = len(cols[0])
            else:
                ncols = 2 + int(r.next() % 4)
                shapes = [gc.SHAPES[int(r.next() % len(gc.SHAPES))] for _ in range(ncols)]
                n = [1, 2, 3, 4, 5, 8, 20, 60, 200][int(r.next() % 9)]
                cols = []
                for j, s in enumerate(shapes):
                    cols.append(gc.column(r, s, n, cols[j - 1] if j else []))
            names = ["c%d_%s" % (j, s) for j, s in enumerate(shapes)]
            rows = [[cols[j][i] for j in range(ncols)] for i in range(n)]
            num = gc.numeric_columns(cols)
            level = LEVELS[int(r.next() % len(LEVELS))]
            q = ot.z_quantile(level)
            for mode in ("pairwise", "listwise"):
                p = run(a.bin, tmp, names, rows, ["--missing", mode, "--level", level])
                if len(num) < 2:
                    doc = json.loads(p.stdout or "{}")
                    if p.returncode != 2 or doc.get("error", {}).get("rule") != "stats.corr-columns":
                        bad.append((fi, "fewer than two numeric columns: want stats.corr-columns, got exit %d %s" % (p.returncode, p.stdout[:100])))
                    continue
                if p.returncode != 0:
                    keep = pathlib.Path(tempfile.gettempdir()) / ("gate_ci_corr_failed_%d.csv" % fi)
                    keep.write_text((tmp / "f.csv").read_text())
                    print("gate_ci_corr: exit %d on file %d (%s): %s (the file is %s)" % (p.returncode, fi, mode, p.stdout[:200], keep))
                    return 2
                doc = json.loads(p.stdout)["data"]
                if doc.get("ci") != "corr" or doc.get("missing") != mode or not gd.same(doc.get("level"), ot.rounded(Fraction(*[ot.level_parts(level)[0], 10 ** ot.level_parts(level)[1]]))) or "pairs" in doc or "adjustment" in doc:
                    bad.append((fi, "the header: %r" % {k: doc.get(k) for k in ("ci", "missing", "level", "rows")}))
                vals = {j: [gc.parse(x) for x in cols[j]] for j in num}
                used_all = [i for i in range(n) if all(vals[j][i] is not None for j in num)]
                want_pairs = [(x, y) for xi, x in enumerate(num) for y in num[xi + 1:]]
                got = doc["intervals"]
                if [(g["x"], g["y"]) for g in got] != [(names[x], names[y]) for x, y in want_pairs]:
                    bad.append((fi, "%s: the pairs are %r" % (mode, [(g["x"], g["y"]) for g in got])))
                    continue
                for g, (x, y) in zip(got, want_pairs):
                    used = used_all if mode == "listwise" else [i for i in range(n) if vals[x][i] is not None and vals[y][i] is not None]
                    w = gc.expected([vals[x][i] for i in used], [vals[y][i] for i in used])
                    pairs += 1
                    tag = "%s %s~%s level %s" % (mode, names[x], names[y], level)
                    if g["n"] != w["n"]:
                        bad.append((fi, "%s n: want %r got %r" % (tag, w["n"], g["n"])))
                    if not gd.same(g.get("estimate"), w["r"]):
                        bad.append((fi, "%s estimate: want %r got %r" % (tag, w["r"], g.get("estimate"))))
                    xs_ = [vals[x][i] for i in used]
                    ys_ = [vals[y][i] for i in used]
                    n_, N_, Lx_, Ly_ = gc.exact_pair(xs_, ys_)
                    if Lx_ == 0 or Ly_ == 0:
                        want_reason = "constant-column"
                    elif N_ * N_ == Lx_ * Ly_:
                        want_reason = "perfect-correlation"
                    elif n_ < 4:
                        want_reason = "too-few-rows"
                    else:
                        want_reason = None
                    if g.get("reason") != want_reason:
                        bad.append((fi, "%s reason: want %r got %r" % (tag, want_reason, g.get("reason"))))
                    if want_reason is not None:
                        if g.get("interval") is not None:
                            bad.append((fi, "%s: want a null interval, got %r" % (tag, g["interval"])))
                        continue
                    iv = g["interval"]
                    if iv is None or iv["of"] != "correlation" or iv["method"] != "fisher-z" or iv["estimate"] != g["estimate"]:
                        bad.append((fi, "%s interval object: %r" % (tag, iv)))
                        continue
                    intervals += 1
                    z, qse, lo, hi = fisher(g["estimate"], w["n"], q)
                    for key, want in (("low", lo), ("high", hi)):
                        if not within(iv[key], want, z, qse):
                            bad.append((fi, "%s %s: want %s got %r" % (tag, key, mpmath.nstr(want, 17), iv[key])))
                        else:
                            scale = max(abs(want), (1 - want * want) * max(abs(z), qse))
                            worst = max(worst, float(abs(mpmath.mpf(iv[key]) - want) / scale / mpmath.mpf("1e-12")))
                    # controls
                    if w["n"] > 4:
                        _, qse2, lo2, hi2 = fisher(g["estimate"], w["n"], q, df_off=2)
                        control("n-minus-2", not (within(iv["low"], lo2, z, qse2) and within(iv["high"], hi2, z, qse2)))
                    control("no-tanh", not (within(float(z - qse), lo, z, qse) and within(float(z + qse), hi, z, qse)))
            if n > 1 and len(num) >= 2:
                p0 = run(a.bin, tmp, names, rows, ["--level", level])
                p2 = run(a.bin, tmp, names, shuffled(r, rows), ["--level", level])
                if p0.stdout != p2.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
        # pairs whose r is within 1e-9 of 1 (z beyond 10) and of -1, and r of the size 1e-9 (the digits of a small r)
        extra = []
        xs = [float(i) for i in range(1, 101)]
        extra.append(("near +1", xs, [x + 1e-4 * (i % 7) for i, x in enumerate(xs)]))
        extra.append(("near -1", xs, [-x - 1e-4 * (i % 5) for i, x in enumerate(xs)]))
        extra.append(("near 0", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], [1.0, 3.0, 2.0, 5.0, 4.0, 7.0, 6.0, 8.0 + 1e-9]))
        extra.append(("small r", xs, [((i * 7919) % 101) + 1e-12 * x for i, x in enumerate(xs)]))
        for k, (what, xs_, ys_) in enumerate(extra):
            rows = [[repr(x_), repr(y_)] for x_, y_ in zip(xs_, ys_)]
            for level in ("0.95", "0.999"):
                p = run(a.bin, tmp, ["a", "b"], rows, ["--level", level])
                g = json.loads(p.stdout)["data"]["intervals"][0]
                n_, N_, Lx_, Ly_ = gc.exact_pair(xs_, ys_)
                w = gc.expected(xs_, ys_)
                pairs += 1
                if not gd.same(g.get("estimate"), w["r"]):
                    bad.append((3000 + k, "%s estimate: want %r got %r" % (what, w["r"], g.get("estimate"))))
                if g.get("interval") is None:
                    bad.append((3000 + k, "%s: no interval (reason %r)" % (what, g.get("reason"))))
                    continue
                intervals += 1
                z, qse, lo, hi = fisher(g["estimate"], n_, ot.z_quantile(level))
                for key, want in (("low", lo), ("high", hi)):
                    if not within(g["interval"][key], want, z, qse):
                        bad.append((3000 + k, "%s level %s %s: want %s got %r" % (what, level, key, mpmath.nstr(want, 17), g["interval"][key])))
        # the other methods are refused
        names = ["a", "b"]
        rows = [["1", "2"], ["2", "1"], ["3", "5"], ["4", "4"], ["5", "7"]]
        for m in ("spearman", "kendall"):
            p = run(a.bin, tmp, names, rows, ["--method", m])
            doc = json.loads(p.stdout) if p.stdout else {}
            if p.returncode != 2 or doc.get("error", {}).get("rule") != "args.bad-value":
                bad.append((-2, "--method %s: want exit 2 args.bad-value, got %d" % (m, p.returncode)))
    for k, (outside, total) in sorted(ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_ci_corr: %d pairs, %d intervals against mpmath (largest error %.3g of the bound), %d problems; controls outside the bound: %s" % (
        pairs, intervals, worst, len(bad), ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(ctl.items()))))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
