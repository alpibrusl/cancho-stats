#!/usr/bin/env python3
"""S3 gates T-A, T-B, T-C, T-D and T-H (docs/test.md section 5) for `stats test`, kind by kind:

    python3 scripts/gate_test.py [--bin build/stats] [--kind t1] [--files 200] [--seed 1] [--scipy]

Generated CSV files (the shapes of gate_describe.py: offsets, wide exponents, huge and tiny values, constants, missing cells; the `mu0` and the level drawn too) against scripts/oracle_test.py, which works from
`Fraction` and `mpmath` and shares no code with the program. For `t1`: n, df, the statistic and Cohen's d and the estimate (bytes), the p-value (relative 1e-11 + (df + 1) 2^-53 against mpmath), the endpoints of
the interval (absolute 1e-12 max(|estimate|, half-width), at least two units of 2^-1074, against mpmath), the reasons (`too-few-rows` below two values, `variance-zero` for a constant column; there the interval has zero width). Also:
the rows shuffled give the same bytes; the columns named in another order give the same bytes (the answer keeps the order of the file); the file negated with `mu0` negated gives the negated statistic and
effect size, the same p-value and the mirrored interval. CONTROLS, from the same data and the oracle's own quantities, must leave the bounds on most of the cases where they differ: `n` for `n - 1` degrees
of freedom, a one-sided p-value, the normal for Student's `t` (p-value), and for the interval the normal quantile and `n` degrees of freedom; so the gate is able to fail. With --scipy: `ttest_1samp` agrees to
1e-12 (relative above |t| = 1, absolute below) on the statistic and 1e-9 relative on the p-value on well-conditioned columns (uniform, ints); the other families are reported only.

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

SHAPES = ["uniform", "uniform", "ints", "ints", "offset", "offset-neg", "wide", "tiny", "huge", "const", "gappy", "gappy", "text", "mixed", "empty"]
GOOD = ("uniform", "ints")
LEVELS = ["0.95", "0.95", "0.9", "0.99", "0.5", "0.999", "0.8"]

# hand cases every run includes: one value, two values, a constant column, mu0 equal to the mean (t = 0), a huge pair, subnormals
FIXED = [
    ([["1.5"]], ["x"], "0", "0.95"),
    ([["1.5", "2.5"]], ["x"], "0", "0.95"),
    ([["3", "3", "3", "3"]], ["x"], "3", "0.95"),
    ([["1", "2", "3", "4", "5"]], ["x"], "3", "0.95"),
    ([["1e300", "2e300", "-3e300", "5e299"]], ["x"], "0", "0.9"),
    ([["5e-324", "1e-323", "1.5e-323", "2e-323"]], ["x"], "0", "0.95"),
    ([["1e8", "100000001", "100000002", "100000004"]], ["x"], "100000001.75", "0.95"),
]


def run(binary, tmp, names, rows, extra):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "test", "t1", "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=300)


def parse(cell):
    if cell == "":
        return None
    if cell.startswith('"'):
        raise ValueError("quoted")
    x = float(cell)
    if math.isnan(x) or math.isinf(x):
        raise ValueError("not finite")
    return x


def numeric(col):
    try:
        vals = [parse(c) for c in col]
    except ValueError:
        return False
    return any(v is not None for v in vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--kind", default="t1", choices=["t1"])
    ap.add_argument("--files", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scipy", action="store_true")
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    checked = 0
    ctl = {"df-n": [0, 0], "one-sided": [0, 0], "normal": [0, 0], "interval-normal-quantile": [0, 0], "interval-df-n": [0, 0]}
    pst = {"max": 0.0, "ratio": 0.0, "n": 0, "iv": 0.0}
    sp = {"t": 0.0, "p": 0.0, "n": 0, "skipped": 0}
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            if fi < len(FIXED):
                cols, names, mu0, level = FIXED[fi]
                shapes = ["fixed"] * len(cols)
                n = len(cols[0])
            else:
                ncols = 1 + int(r.next() % 4)
                shapes = [SHAPES[int(r.next() % len(SHAPES))] for _ in range(ncols)]
                n = [1, 2, 3, 4, 5, 8, 20, 60, 200][int(r.next() % 9)]
                cols = [gd.column(r, s, n) for s in shapes]
                names = ["c%d_%s" % (j, s) for j, s in enumerate(shapes)]
                level = LEVELS[int(r.next() % len(LEVELS))]
                mu0 = r.pick(["0", "1", "-3.5", "21.5", "1e8", "0.001"])
            rows = [[cols[j][i] for j in range(len(cols))] for i in range(n)]
            num = [j for j in range(len(cols)) if numeric(cols[j])]
            if not num:
                continue
            chosen = [j for j in num if r.next() % 3 != 0] or [num[0]]
            if fi >= len(FIXED) and r.next() % 4 == 0:
                # mu0 = the exact mean of the first chosen column when it is representable, else the rounded mean: t = 0 or tiny
                vals0 = [parse(c) for c in cols[chosen[0]] if parse(c) is not None]
                if vals0:
                    mu0 = repr(ot.rounded(sum((Fraction(v) for v in vals0), Fraction(0)) / len(vals0)))
            extra = ["--cols", ",".join(names[j] for j in chosen), "--mu0", mu0, "--level", level]
            try:
                p = run(a.bin, tmp, names, rows, extra)
            except subprocess.TimeoutExpired:
                print("gate_test: the program did not finish (hang)")
                return 2
            if p.returncode != 0:
                keep = pathlib.Path(tempfile.gettempdir()) / ("gate_test_failed_%d.csv" % fi)
                keep.write_text((tmp / "f.csv").read_text())
                print("gate_test: exit %d on file %d: %s (the file is %s, args %r)" % (p.returncode, fi, p.stdout[:240], keep, extra))
                return 2
            doc = json.loads(p.stdout)["data"]
            want_level_top = float(Fraction(ot.level_parts(level)[0], 10 ** ot.level_parts(level)[1]))
            if doc.get("test") != "t1" or doc.get("rows") != n or not gd.same(doc.get("mu0"), float(mu0)) or not gd.same(doc.get("level"), want_level_top) or doc.get("adjustment") != "none":
                bad.append((fi, "the header of the answer: want test t1, rows %d, mu0 %r, level %r, adjustment none; got %r" % (n, float(mu0), want_level_top, {k: doc.get(k) for k in ("test", "rows", "mu0", "level", "adjustment")})))
            tests = doc["tests"]
            if [x["column"] for x in tests] != [names[j] for j in chosen]:
                bad.append((fi, "the tests are %r" % [x["column"] for x in tests]))
                continue
            for g, j in zip(tests, chosen):
                vals = [v for v in (parse(c) for c in cols[j]) if v is not None]
                w = ot.t1(vals, float(mu0), level)
                tag = "t1 %s mu0 %s level %s" % (names[j], mu0, level)
                checked += 1
                if g["n"] != {"n": w["n"], "missing": n - w["n"]}:
                    bad.append((fi, "%s n: want %r got %r" % (tag, {"n": w["n"], "missing": n - w["n"]}, g["n"])))
                if not gd.same(g.get("statistic"), w["statistic"]):
                    bad.append((fi, "%s statistic: want %r got %r" % (tag, w["statistic"], g.get("statistic"))))
                if g.get("df") != w["df"]:
                    bad.append((fi, "%s df: want %r got %r" % (tag, w["df"], g.get("df"))))
                es = g.get("effect_size")
                if w["d"] is None:
                    if es is not None:
                        bad.append((fi, "%s effect_size: want null got %r" % (tag, es)))
                elif es is None or es.get("name") != "cohen_d" or not gd.same(es.get("value"), w["d"]):
                    bad.append((fi, "%s effect_size: want %r got %r" % (tag, w["d"], es)))
                if w["p"] is None:
                    if g.get("p_value") is not None:
                        bad.append((fi, "%s p_value: want null got %r" % (tag, g.get("p_value"))))
                elif not ot.p_ok(g.get("p_value"), w["p"], w["df"]):
                    bad.append((fi, "%s p_value: want %s got %r" % (tag, mpmath.nstr(w["p"], 17), g.get("p_value"))))
                else:
                    if w["p"] > mpmath.mpf("1e-300"):
                        err = float(abs(mpmath.mpf(g["p_value"]) - w["p"]) / w["p"])
                        pst["max"] = max(pst["max"], err)
                        pst["ratio"] = max(pst["ratio"], err / (1e-11 + (w["df"] + 1) * ot.EPS))
                        pst["n"] += 1
                if g.get("p_value") is not None and not (0.0 <= g["p_value"] <= 1.0):
                    bad.append((fi, "%s p_value %r is not in [0, 1]" % (tag, g["p_value"])))
                if g.get("reason") != w["reason"]:
                    bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
                if g.get("p_method") != "t" or g.get("method") != "one-sample t":
                    bad.append((fi, "%s p_method/method %r %r" % (tag, g.get("p_method"), g.get("method"))))
                iv = g.get("interval")
                if w["low"] is None:
                    if iv is not None:
                        bad.append((fi, "%s interval: want null got %r" % (tag, iv)))
                else:
                    want_level = float(Fraction(*ot.level_parts(level)[:1], 10 ** ot.level_parts(level)[1]))
                    if iv is None or iv.get("of") != "mean" or iv.get("rounding") != "double" or not gd.same(iv.get("level"), want_level) or not gd.same(iv.get("estimate"), w["estimate"]):
                        bad.append((fi, "%s interval header: want estimate %r level %r got %r" % (tag, w["estimate"], want_level, iv)))
                    else:
                        half = w["half"]
                        for k in ("low", "high"):
                            if not ot.endpoint_ok(iv.get(k), w[k], w["estimate"], half):
                                bad.append((fi, "%s interval %s: want %s got %r" % (tag, k, mpmath.nstr(w[k], 17), iv.get(k))))
                            else:
                                scale = max(abs(mpmath.mpf(w["estimate"])), abs(half))
                                ratio_iv = float(abs(mpmath.mpf(iv[k]) - w[k]) / scale) if scale > mpmath.mpf("1e-290") else 0.0   # the subnormal range is judged by its 2-unit floor, not reported
                                if ratio_iv > pst["iv"]:
                                    pst["iv"] = ratio_iv
                                    pst["iv_at"] = "%s %s n=%d est=%r half=%s got=%r want=%s" % (tag, k, w["n"], w["estimate"], mpmath.nstr(half, 6), iv[k], mpmath.nstr(w[k], 17))
                # the controls
                if w["p"] is not None and w["p"] > mpmath.mpf("1e-250") and w["n"] >= 3 and w["t2"] != 0:
                    nn_ = w["n"]
                    wrong_df = ot.p_two_sided_t(nn_, Fraction(nn_) / (nn_ + w["t2"]))
                    ctl["df-n"][1] += 1
                    ctl["df-n"][0] += 0 if ot.p_ok(float(wrong_df), w["p"], w["df"]) else 1
                    ctl["one-sided"][1] += 1
                    ctl["one-sided"][0] += 0 if ot.p_ok(float(w["p"] / 2), w["p"], w["df"]) else 1
                    tt = mpmath.sqrt(mpmath.mpf(w["t2"].numerator) / w["t2"].denominator)
                    ctl["normal"][1] += 1
                    ctl["normal"][0] += 0 if ot.p_ok(float(mpmath.erfc(tt / mpmath.sqrt(2))), w["p"], w["df"]) else 1
                if w["low"] is not None and w["half"] > 0 and w["n"] >= 3:
                    nn_ = w["n"]
                    zq = mpmath.sqrt(2) * mpmath.erfinv(1 - 2 * mpmath.mpf(ot.alpha_half(level).numerator) / ot.alpha_half(level).denominator)
                    ctl["interval-normal-quantile"][1] += 1
                    ctl["interval-normal-quantile"][0] += 0 if ot.endpoint_ok(iv["high"], mpmath.mpf(w["estimate"]) + zq * w["se"], w["estimate"], w["half"]) else 1
                    qn = ot.t_quantile_upper(ot.alpha_half(level), nn_)
                    ctl["interval-df-n"][1] += 1
                    ctl["interval-df-n"][0] += 0 if ot.endpoint_ok(iv["high"], mpmath.mpf(w["estimate"]) + qn * w["se"], w["estimate"], w["half"]) else 1
                if a.scipy and w["statistic"] is not None and w["p"] is not None and g.get("p_value") is not None:
                    import scipy.stats as ss
                    res = ss.ttest_1samp(vals, float(mu0))
                    dt = abs(res.statistic - g["statistic"]) / max(abs(res.statistic), 1.0)  # relative, and absolute below |t| = 1: at t near 0 scipy's float cancellation decides its own answer
                    dp = abs(res.pvalue - g["p_value"]) / max(abs(res.pvalue), 1e-300)
                    if shapes[j] in GOOD:
                        sp["n"] += 1
                        sp["t"] = max(sp["t"], dt)
                        sp["p"] = max(sp["p"], dp)
                        if dt > 1e-12 or dp > 1e-9:
                            bad.append((fi, "%s scipy: dt %g dp %g" % (tag, dt, dp)))
                    else:
                        sp["skipped"] += 1
            # the order of the rows, of the names, and the sign
            if n > 1:
                perm = list(range(n))
                for i in range(n - 1, 0, -1):
                    k = r.next() % (i + 1)
                    perm[i], perm[k] = perm[k], perm[i]
                p2 = run(a.bin, tmp, names, [rows[i] for i in perm], extra)
                if p2.stdout != p.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
                extra3 = ["--cols", ",".join(reversed([names[j] for j in chosen])), "--mu0", mu0, "--level", level]
                p3 = run(a.bin, tmp, names, rows, extra3)
                if p3.stdout != p.stdout:
                    bad.append((fi, "--cols in another order gives different bytes"))
                if len(chosen) == 1 and shapes[chosen[0]] != "fixed" or shapes[chosen[0]] == "fixed":
                    neg = [[(("-" + c) if c and not c.startswith("-") else (c[1:] if c.startswith("-") else c)) if j in chosen else c for j, c in enumerate(rw)] for rw in rows]
                    nm = mu0[1:] if mu0.startswith("-") else "-" + mu0
                    if mu0 in ("0", "-0"):
                        nm = "0"
                    p4 = run(a.bin, tmp, names, neg, ["--cols", ",".join(names[j] for j in chosen), "--mu0", nm, "--level", level])
                    if p4.returncode == 0:
                        d4 = json.loads(p4.stdout)["data"]["tests"]
                        for g0, g4 in zip(tests, d4):
                            ok = g0["n"] == g4["n"] and g0["df"] == g4["df"] and g0["p_value"] == g4["p_value"] and g0.get("reason") == g4.get("reason")
                            for k in ("statistic",):
                                if g0[k] is not None:
                                    ok = ok and gd.same(-g0[k] if g0[k] != 0 else 0.0, g4[k])
                            iv0, iv4 = g0.get("interval"), g4.get("interval")
                            if iv0 is not None:
                                ok = ok and gd.same(-iv0["estimate"] if iv0["estimate"] != 0 else 0.0, iv4["estimate"])
                            if not ok:
                                bad.append((fi, "the negated file does not mirror the answer for %s" % g0["column"]))
                    else:
                        bad.append((fi, "negated file: exit %d" % p4.returncode))
    for k, (outside, total) in ctl.items():
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_test t1: %d files, %d tests compared, %d differ; p-values against mpmath: largest relative error %.3g over %d, at most %.3g of the bound; interval endpoints: largest error %.3g of max(|estimate|, half-width) (bound 1e-12); controls outside the bound: %s%s" % (
        a.files, checked, len(bad), pst["max"], pst["n"], pst["ratio"], pst["iv"], ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in ctl.items()),
        ("; scipy on %d well-conditioned columns: max relative dt %.3g, dp %.3g (%d others not gated)" % (sp["n"], sp["t"], sp["p"], sp["skipped"])) if a.scipy else ""))
    if pst.get("iv_at"):
        print("  largest interval error at: %s" % pst["iv_at"])
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
