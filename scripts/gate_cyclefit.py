#!/usr/bin/env python3
"""S6.2 gates G-F1 to G-F4 for `stats cycle fit` (docs/estimation.md sections 4 and 7):

    python3 scripts/gate_cyclefit.py [--bin build/stats] [--cases 150] [--seed 1]

An independent implementation recomputes every series in exact rational arithmetic: the weighted design is built in doubles (square root of the weight, the polynomial of the scaled time, the sine and cosine of 2 pi t / period),
and the least squares of THAT design is solved exactly (normal equations in fractions): the coefficients, the residual sum of squares, the 2 x 2 block of the inverse for the cycle pair, the Wald statistic F with its 1e-18
ridge, and the shrinkage max(0, 1 - 2 / F). The program's numbers must equal them to the RELATIVE BOUND of docs/estimation.md section 7 (1e-9, with a floor of 1e-9 of the largest coefficient); the maximum gap is printed. Also:
the rows in another order give the same bytes (the days are sorted), a series of fewer than 3 (order + 3) days is `too-short` with zeros and exactly that many is fitted, `--shrink none` gives the raw fit, the rank-deficient
design (a period of one day on whole days) is `degenerate`, and every refusal is reached with its tag. CONTROLS: the shrinkage with the sign of 2 / F changed, the polynomial in absolute time and the divisor n instead of
n - p must each leave the answer on at least a tenth of the cases.

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

BOUND = 1e-9


def run(binary, tmp, text, argv):
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "cycle", "fit"] + list(argv) + ["--root", str(tmp), "f.csv"], capture_output=True, text=True, timeout=600)


def solve(A, b):
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None:
            return None
        M[c], M[piv] = M[piv], M[c]
        for r in range(n):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [M[r][k] - f * M[c][k] for k in range(n + 1)]
    return [M[i][n] / M[i][i] for i in range(n)]


def inverse(A):
    n = len(A)
    cols = []
    for j in range(n):
        e = [Fraction(int(i == j)) for i in range(n)]
        cols.append(solve(A, e))
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def exact_fit(days, ys, ws, period, order, shrink="james-stein", wrong=None):
    """(status, a_raw, b_raw, f, shrink, a, b) of one series as exact fractions, or the status alone"""
    if wrong == "order":
        order = order + 1 if order < 3 else order - 1
    p = order + 3
    pairs = sorted((d, y, w) for d, y, w in zip(days, ys, ws) if w > 0)
    n = len(pairs)
    if n < 3 * (order + 3):
        return ("too-short",)
    t0, t1 = pairs[0][0], pairs[-1][0]
    span = max(float(t1 - t0), 1.0)
    X, Y = [], []
    for d, y, w in pairs:
        sw = math.sqrt(w)
        own = (float(d) - float(t0)) / span
        row = [(own ** k) * sw for k in range(order + 1)]
        ang = 6.283185307179586 * d / period
        row += [math.sin(ang) * sw, math.cos(ang) * sw]
        X.append([Fraction(v) for v in row])
        Y.append(Fraction(y * sw))
    A = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(p)] for a in range(p)]
    rhs = [sum(X[i][a] * Y[i] for i in range(n)) for a in range(p)]
    c = solve(A, rhs)
    if c is None:
        return ("degenerate",)
    rss = sum((Y[i] - sum(X[i][a] * c[a] for a in range(p))) ** 2 for i in range(n))
    dof = max(n - p, 1) if wrong != "n" else n
    s2 = rss / dof
    inv = inverse(A)
    ar, br = c[p - 2], c[p - 1]
    ridge = Fraction("1e-18")
    m00, m01, m11 = s2 * inv[p - 2][p - 2] + ridge, s2 * inv[p - 2][p - 1], s2 * inv[p - 1][p - 1] + ridge
    det = m00 * m11 - m01 * m01
    f = (ar * (m11 * ar - m01 * br) + br * (m00 * br - m01 * ar)) / det if det != 0 else Fraction(0)
    if shrink == "none":
        sh = Fraction(1)
    elif f > 0:
        sh = max(Fraction(0), 1 - 2 / f) if wrong != "sign" else max(Fraction(0), 1 + 2 / f)
    else:
        sh = Fraction(0)
    return ("fitted", ar, br, f, sh, sh * ar, sh * br)


def close(got, want, scale):
    w = float(want)
    return abs(got - w) <= BOUND * (abs(w) + 1e-9 * scale)


def gap(got, want, scale):
    w = float(want)
    return abs(got - w) / (abs(w) + 1e-9 * scale)


def parse(r):
    d = json.loads(r.stdout)
    return d["ok"], d.get("data"), d.get("error")


def make_file(rng, nseries):
    """text of a long-format file and the series as {key: (days, ys, ws, has_weight)}"""
    series = {}
    order_rows = []
    use_w = rng.next() % 2 == 0
    for k in range(nseries):
        key = rng.pick(["north", "south", "FR", "es-01", "k%d" % k, "Z"]) + str(k)
        n = rng.pick([14, 15, 16, 20, 40, 90, 200, 366])
        gaps = rng.next() % 3 == 0
        days, d = [], rng.pick([0, 0, 7, 100])
        for i in range(n):
            days.append(d)
            d += 1 + (rng.next() % 3 if gaps else 0)
        period = None
        amp = rng.pick([0.0, 0.1, 0.5, 2.0])
        ph = rng.u() * 6
        ys = [round(amp * math.sin(2 * math.pi * x / 60 + ph) + 0.01 * x + (rng.u() - 0.5) * rng.pick([0.1, 1.0, 5.0]), 6) for x in days]
        ws = [round(0.2 + rng.u() * 3, 4) if use_w else 1.0 for _ in days]
        if use_w and rng.next() % 6 == 0:
            ws[rng.next() % n] = 0.0                      # a zero weight: left out
        series[key] = (days, ys, ws)
        for i in range(n):
            order_rows.append((key, days[i], ys[i], ws[i]))
    return series, order_rows, use_w


def text_of(rows, use_w):
    head = "k,day,y" + (",w" if use_w else "")
    return head + "\n" + "\n".join(",".join(str(v) for v in (r if use_w else r[:3])) for r in rows) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--cases", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = SM64(a.seed)
    bad = []
    worst = {"a_raw": 0.0, "b_raw": 0.0, "f": 0.0, "shrink": 0.0, "a": 0.0, "b": 0.0}
    differs = {"sign": 0, "order": 0, "n": 0}
    degenerate = 0
    done = 0
    nser = 0
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        for k in range(a.cases):
            series, rows, use_w = make_file(rng, rng.pick([1, 2, 3, 5]))
            order = rng.pick([0, 1, 2, 2, 2, 3])
            period = rng.pick([21.0, 30.5, 45.0, 60.0, 90.0, 91.25, 180.0, 365.0, 7.3])
            shrink = rng.pick(["james-stein", "james-stein", "none"])
            argv = ["--day", "day", "--y", "y", "--series", "k", "--period", repr(period), "--order", str(order), "--shrink", shrink] + (["--weight", "w"] if use_w else [])
            r = run(a.bin, tmp, text_of(rows, use_w), argv)
            ok, data, err = parse(r)
            if not ok:
                bad.append(("case %d refused" % k, err))
                continue
            # the rows reversed: the same bytes
            r2 = run(a.bin, tmp, text_of(list(reversed(rows)), use_w), argv)
            if r2.stdout != r.stdout:
                bad.append(("case %d order" % k,))
            got = {s["series"]: s for s in data["series"]}
            if set(got) != set(series):
                bad.append(("case %d series" % k, sorted(got), sorted(series)))
                continue
            for key, (days, ys, ws) in series.items():
                want = exact_fit(days, ys, ws, period, order, shrink)
                g = got[key]
                nser += 1
                if want[0] != g["status"] and not (want[0] == "fitted" and g["status"] == "degenerate"):
                    bad.append(("case %d %s status" % (k, key), g["status"], want[0]))
                    continue
                if want[0] == "too-short":
                    if (g["a"], g["b"], g["a_raw"], g["b_raw"], g["f"], g["shrink"]) != (0.0,) * 6:
                        bad.append(("case %d %s too-short not zero" % (k, key), g))
                    continue
                if g["status"] == "degenerate":
                    degenerate += 1
                    continue
                if want[0] == "degenerate":
                    continue
                scale = max(abs(float(want[1])), abs(float(want[2])), 1e-9)
                for name, w_, got_ in (("a_raw", want[1], g["a_raw"]), ("b_raw", want[2], g["b_raw"]), ("f", want[3], g["f"]), ("shrink", want[4], g["shrink"]), ("a", want[5], g["a"]), ("b", want[6], g["b"])):
                    sc = scale if name != "f" else max(abs(float(want[3])), 1e-9)
                    worst[name] = max(worst[name], gap(got_, w_, sc))
                    if not close(got_, w_, sc):
                        bad.append(("case %d %s %s" % (k, key, name), got_, float(w_)))
                if g["days"] != sum(1 for w in ws if w > 0):
                    bad.append(("case %d %s days" % (k, key), g["days"]))
                if shrink == "james-stein":
                    for name in differs:
                        w2 = exact_fit(days, ys, ws, period, order, shrink, wrong=name)
                        if w2[0] == "fitted" and not (close(g["a"], w2[5], scale) and close(g["b"], w2[6], scale) and close(g["shrink"], w2[4], max(abs(float(w2[4])), 1e-9))):
                            differs[name] += 1
            done += 1
        if nser and degenerate * 5 > nser:
            bad.append(("%d of %d series reported degenerate: the refusal is hiding errors" % (degenerate, nser),))
        for name, count in differs.items():
            if nser and count * 10 < nser // 2:
                bad.append(("control %s left the answer on %d of about %d series" % (name, count, nser // 2),))
        # ---- G-F4: the length of a too-short series, exactly at the edge
        for order in (0, 2, 3):
            need = 3 * (order + 3)
            for n, status in ((need - 1, "too-short"), (need, "fitted")):
                rowsx = [("s", i, round(math.sin(i / 3.0) + 0.1 * (i % 5), 6), 1.0) for i in range(n)]
                r = run(a.bin, tmp, text_of(rowsx, False), ["--day", "day", "--y", "y", "--series", "k", "--period", "7.3", "--order", str(order)])
                ok, data, err = parse(r)
                if not ok or data["series"][0]["status"] != status:
                    bad.append(("edge order %d n %d" % (order, n), ok and data["series"][0]["status"], err))
        # the period of one day on whole days: the sine is zero and the cosine one, so the design is rank deficient
        rowsx = [("s", i, round(math.sin(i / 3.0), 6), 1.0) for i in range(40)]
        r = run(a.bin, tmp, text_of(rowsx, False), ["--day", "day", "--y", "y", "--series", "k", "--period", "1"])
        ok, data, err = parse(r)
        if not ok or data["series"][0]["status"] != "degenerate":
            bad.append(("period 1", data and data["series"][0]["status"], err))
        # a series key that is missing (an empty cell) is counted and left out
        r = run(a.bin, tmp, "k,day,y\n,1,2\n" + "\n".join("s,%d,%r" % (i, math.sin(i)) for i in range(20)) + "\n", ["--day", "day", "--y", "y", "--series", "k", "--period", "7"])
        ok, data, err = parse(r)
        if not ok or data["series_missing"] != 1 or len(data["series"]) != 1:
            bad.append(("series_missing", data, err))
        # ---- refusals
        base = "k,day,y\n" + "\n".join("s,%d,%r" % (i, math.sin(i)) for i in range(20)) + "\n"
        argv0 = ["--day", "day", "--y", "y", "--series", "k", "--period", "7"]
        refusals = [
            ("stats.test-columns", base, ["--day", "day", "--y", "y", "--series", "k"]),
            ("args.bad-value", base, ["--day", "day", "--y", "y", "--series", "k", "--period", "abc"]),
            ("args.bad-value", base, ["--day", "day", "--y", "y", "--series", "k", "--period", "0"]),
            ("args.bad-value", base, argv0 + ["--order", "4"]),
            ("stats.no-such-column", base, ["--day", "nope", "--y", "y", "--series", "k", "--period", "7"]),
            ("stats.non-numeric", base, ["--day", "k", "--y", "y", "--series", "k", "--period", "7"]),
            ("stats.cycle.duplicate-day", base + "s,3,1\n", argv0),
            ("stats.cycle.day-not-integer", "k,day,y\ns,1.5,2\ns,2,3\n", argv0),
            ("stats.cycle.day-not-integer", "k,day,y\ns,-1,2\ns,2,3\n", argv0),
            ("stats.cycle.weight-negative", "k,day,y,w\ns,1,2,-1\ns,2,3,1\n", argv0 + ["--weight", "w"]),
        ]
        for tag, text, av in refusals:
            r = run(a.bin, tmp, text, av)
            try:
                ok, data, err = parse(r)
            except ValueError:
                bad.append(("refusal %s: no JSON" % tag, r.stdout[:100], r.stderr[:100]))
                continue
            if ok or err["rule"] != tag:
                bad.append(("refusal %s" % tag, ok, err and err["rule"]))
    print("cycle fit: %d files, %d series; worst relative gap a_raw %.2e b_raw %.2e f %.2e shrink %.2e a %.2e b %.2e; controls differed on %s series; %d reported degenerate" % (
        done, nser, worst["a_raw"], worst["b_raw"], worst["f"], worst["shrink"], worst["a"], worst["b"], differs, degenerate))
    if bad:
        for b in bad[:20]:
            print("FAIL", *b)
        print("%d differences" % len(bad))
        return 1
    print("gate_cyclefit: ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
