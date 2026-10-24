#!/usr/bin/env python3
"""S3.4 gates T-A, T-B, T-C and T-H (docs/test.md sections 5 and 11) for `stats test prop`:

    python3 scripts/gate_test4.py [--bin build/stats] [--files 200] [--seed 1]

Generated CSV files (an outcome column with a success label among others, a group column with two labels, empty cells in both; trials from 1 to 3,000; proportions of 0, 1 and everything between; hand cases x = 0, x = n,
one trial, a proportion equal to p0, p0 of nine decimals) against scripts/oracle_test.py, which works from `Fraction` and `mpmath` and shares no code with the program. One sample (`--success S --p0 P`): the statistic and
the difference p - p0 BIT FOR BIT (each one rounding of an exact rational), the successes and trials and the missing rows, the p-value (relative 1e-11 + (z^2 + 1) 2^-53 against mpmath at the exact z), the Wilson
interval (absolute 1e-11 against mpmath at 60 digits), n p0 and n (1 - p0). Two samples (`--group G`): the same for the pooled z and the difference, Newcombe's interval (absolute 1e-11), Cohen's h (absolute 1e-12), the
pooled proportion, the reasons (`variance-zero` when every row or no row is a success). Also: the rows shuffled give the same bytes; the columns named in another order give the same entries in that order; the groups swapped
(the labels renamed so that the byte order flips) negate the statistic, the difference and h and mirror the interval. CONTROLS must leave the bound on a tenth of the cases: Wald's interval, Yates's continuity
correction, the one-sided p-value, Student's t in place of the normal (for the p-value), and for two samples the unpooled z.

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

LEVELS = ["0.95", "0.95", "0.9", "0.99", "0.5", "0.999", "0.8"]
P0S = ["0.5", "0.5", "0.1", "0.25", "0.9", "0.333333333", "0.000000001", "0.999999999", "0.05"]
OUTCOMES = [("yes", ["no"]), ("1", ["0"]), ("a b", ["c", "d"]), ("x,y", ["z"]), ("T", ["F", "U"])]
IV = mpmath.mpf("1e-11")
H = mpmath.mpf("1e-12")


def run(binary, tmp, names, rows, extra):
    text = ",".join(names) + "\n" + "\n".join(",".join(cell(c) for c in r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "test", "prop", "--root", str(tmp)] + list(extra) + ["f.csv"], capture_output=True, text=True, timeout=300)


def cell(v):
    if v is None:
        return ""
    if "," in v or '"' in v:
        return '"' + v.replace('"', '""') + '"'
    return v


class Stat:
    def __init__(self):
        self.ctl = {}
        self.p = {"max": 0.0, "ratio": 0.0, "n": 0}
        self.iv = 0.0
        self.h = 0.0
        self.checked = 0

    def control(self, name, outside):
        c = self.ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0


def near(got, want, tol):
    return got is not None and abs(mpmath.mpf(got) - want) <= tol


def check_p(bad, fi, tag, g, w, st):
    if w["p"] is None:
        if g.get("p_value") is not None:
            bad.append((fi, "%s p_value: want null got %r" % (tag, g.get("p_value"))))
        return
    z2 = float(w["z2"])
    if not ot.p_ok(g.get("p_value"), w["p"], z2):
        bad.append((fi, "%s p_value: want %s got %r" % (tag, mpmath.nstr(w["p"], 17), g.get("p_value"))))
    elif w["p"] > mpmath.mpf("1e-300"):
        err = float(abs(mpmath.mpf(g["p_value"]) - w["p"]) / w["p"])
        st.p["max"] = max(st.p["max"], err)
        st.p["ratio"] = max(st.p["ratio"], err / (1e-11 + (z2 + 1) * ot.EPS))
        st.p["n"] += 1


def check_one(bad, fi, tag, g, w, rows_total, p0, d, level, st):
    st.checked += 1
    if g.get("method") != "one-sample z" or g.get("p_method") != "normal" or g.get("df") is not None:
        bad.append((fi, "%s method/p_method/df: %r" % (tag, {k: g.get(k) for k in ("method", "p_method", "df")})))
    if g["n"] != {"n": w["n"], "successes": w["x"], "missing": rows_total - w["n"]}:
        bad.append((fi, "%s n: want %r got %r" % (tag, {"n": w["n"], "successes": w["x"], "missing": rows_total - w["n"]}, g["n"])))
    if g.get("reason") != w["reason"]:
        bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
    if w["statistic"] is None:
        for k in ("statistic", "p_value", "effect_size", "interval"):
            if g.get(k) is not None:
                bad.append((fi, "%s %s: want null got %r" % (tag, k, g.get(k))))
        return
    if not gd.same(g.get("statistic"), w["statistic"]):
        bad.append((fi, "%s statistic: want %r got %r" % (tag, w["statistic"], g.get("statistic"))))
    es = g.get("effect_size")
    if es is None or es.get("name") != "difference" or es.get("rounding") != "single" or not gd.same(es.get("value"), w["diff"]):
        bad.append((fi, "%s effect_size: want difference %r got %r" % (tag, w["diff"], es)))
    check_p(bad, fi, tag, g, w, st)
    iv = g.get("interval")
    if iv is None or iv.get("of") != "proportion" or iv.get("method") != "wilson" or iv.get("rounding") != "double" or not gd.same(iv.get("estimate"), w["estimate"]):
        bad.append((fi, "%s interval header: %r" % (tag, iv)))
    else:
        for k, want in (("low", w["low"]), ("high", w["high"])):
            if not near(iv.get(k), want, IV):
                bad.append((fi, "%s interval %s: want %s got %r" % (tag, k, mpmath.nstr(want, 17), iv.get(k))))
            else:
                st.iv = max(st.iv, float(abs(mpmath.mpf(iv[k]) - want)))
    a = g.get("assumptions") or {}
    if not gd.same(a.get("n_p0"), w["n_p0"]) and not near(a.get("n_p0"), mpmath.mpf(w["n_p0"]), mpmath.mpf(w["n_p0"]) * mpmath.mpf("2.3e-16")):
        bad.append((fi, "%s n_p0: want %r got %r" % (tag, w["n_p0"], a.get("n_p0"))))


def controls_one(st, g, w, p0, level):
    if w["statistic"] is None or w["p"] is None or not w["p"] > mpmath.mpf("1e-250") or w["n"] < 3:
        return
    z2 = float(w["z2"])
    n, x = w["n"], w["x"]
    st.control("one-sided", not ot.p_ok(float(w["p"] / 2), w["p"], z2))
    sd = mpmath.sqrt(frac_to_mp(n * p0 * (1 - p0)))
    zc = max(abs(frac_to_mp(x - n * p0)) - mpmath.mpf(1) / 2, mpmath.mpf(0)) / sd
    st.control("yates", not ot.p_ok(float(mpmath.erfc(zc / mpmath.sqrt(2))), w["p"], z2))
    z = mpmath.sqrt(frac_to_mp(w["z2"]))
    import scipy.stats as ss
    st.control("student-t", not ot.p_ok(float(2 * ss.t.sf(float(z), n - 1)), w["p"], z2))
    # Wald: p -+ q sqrt(p (1 - p) / n)
    ph = mpmath.mpf(x) / n
    half = w["q"] * mpmath.sqrt(ph * (1 - ph) / n)
    st.control("wald", not (near(g["interval"]["low"], ph - half, IV) and near(g["interval"]["high"], ph + half, IV)))


def frac_to_mp(f):
    return mpmath.mpf(f.numerator) / f.denominator


def check_two(bad, fi, tag, g, w, rows_total, n1, x1, n2, x2, la, lb, st):
    st.checked += 1
    if g.get("method") != "two-sample z (pooled)" or g.get("p_method") != "normal" or g.get("df") is not None:
        bad.append((fi, "%s method/p_method/df: %r" % (tag, {k: g.get(k) for k in ("method", "p_method", "df")})))
    if g["n"] != {"a": {"n": n1, "successes": x1}, "b": {"n": n2, "successes": x2}, "missing": rows_total - n1 - n2}:
        bad.append((fi, "%s n: got %r" % (tag, g["n"])))
    if g.get("levels") != {"a": la, "b": lb}:
        bad.append((fi, "%s levels: want %r got %r" % (tag, {"a": la, "b": lb}, g.get("levels"))))
    if g.get("reason") != w["reason"]:
        bad.append((fi, "%s reason: want %r got %r" % (tag, w["reason"], g.get("reason"))))
    if w["reason"] == "too-few-rows":
        return
    es = g.get("effect_size")
    if es is None or es.get("name") != "cohen_h" or es.get("rounding") != "double" or not near(es.get("value"), w["h"], H):
        bad.append((fi, "%s effect_size: want cohen_h %s got %r" % (tag, mpmath.nstr(w["h"], 17), es)))
    else:
        st.h = max(st.h, float(abs(mpmath.mpf(es["value"]) - w["h"])))
    iv = g.get("interval")
    if iv is None or iv.get("of") != "difference of proportions" or iv.get("method") != "newcombe" or iv.get("rounding") != "double" or not gd.same(iv.get("estimate"), w["diff"]):
        bad.append((fi, "%s interval header: want estimate %r got %r" % (tag, w["diff"], iv)))
    else:
        for k in ("low", "high"):
            if not near(iv.get(k), w[k], IV):
                bad.append((fi, "%s interval %s: want %s got %r" % (tag, k, mpmath.nstr(w[k], 17), iv.get(k))))
            else:
                st.iv = max(st.iv, float(abs(mpmath.mpf(iv[k]) - w[k])))
    if w["statistic"] is None:
        for k in ("statistic", "p_value"):
            if g.get(k) is not None:
                bad.append((fi, "%s %s: want null got %r" % (tag, k, g.get(k))))
        return
    if not gd.same(g.get("statistic"), w["statistic"]):
        bad.append((fi, "%s statistic: want %r got %r" % (tag, w["statistic"], g.get("statistic"))))
    check_p(bad, fi, tag, g, w, st)
    if not gd.same((g.get("assumptions") or {}).get("pooled"), w["pooled"]):
        bad.append((fi, "%s pooled: want %r got %r" % (tag, w["pooled"], (g.get("assumptions") or {}).get("pooled"))))


def controls_two(st, g, w, n1, x1, n2, x2):
    if w["statistic"] is None or w["p"] is None or not w["p"] > mpmath.mpf("1e-250") or n1 < 3 or n2 < 3:
        return
    z2 = float(w["z2"])
    st.control("one-sided", not ot.p_ok(float(w["p"] / 2), w["p"], z2))
    p1, p2 = Fraction(x1, n1), Fraction(x2, n2)
    var = p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2
    if var != 0:
        zu = abs(frac_to_mp(p1 - p2)) / mpmath.sqrt(frac_to_mp(var))
        st.control("unpooled", not ot.p_ok(float(mpmath.erfc(zu / mpmath.sqrt(2))), w["p"], z2))
    nn = n1 + n2
    pp = Fraction(x1 + x2, nn)
    sd = mpmath.sqrt(frac_to_mp(pp * (1 - pp) * (Fraction(1, n1) + Fraction(1, n2))))
    corr = Fraction(1, 2) * (Fraction(1, n1) + Fraction(1, n2))
    zc = max(abs(frac_to_mp(p1 - p2)) - frac_to_mp(corr), mpmath.mpf(0)) / sd
    st.control("yates", not ot.p_ok(float(mpmath.erfc(zc / mpmath.sqrt(2))), w["p"], z2))
    ph = mpmath.sqrt(frac_to_mp(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2))
    d = frac_to_mp(p1 - p2)
    st.control("wald", not (near(g["interval"]["low"], d - w["q"] * ph, IV) and near(g["interval"]["high"], d + w["q"] * ph, IV)))


def shuffled(r, seq):
    perm = list(range(len(seq)))
    for i in range(len(perm) - 1, 0, -1):
        k = r.next() % (i + 1)
        perm[i], perm[k] = perm[k], perm[i]
    return [seq[i] for i in perm]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    st = Stat()
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            succ, others = OUTCOMES[int(r.next() % len(OUTCOMES))]
            n = [1, 2, 3, 5, 10, 30, 100, 400, 3000][int(r.next() % 9)]
            ps = [0.0, 1.0, 0.5, 0.1, 0.9, 0.3][int(r.next() % 6)] if r.next() % 4 == 0 else (r.next() % 1000) / 1000.0
            labels = (succ, others)
            la, lb = ("g1", "g2")
            rows = []
            for _ in range(n):
                u = (r.next() % 100000) / 100000.0
                y = succ if u < ps else others[int(r.next() % len(others))]
                # the second group has another proportion
                grp = la if r.next() % 2 else lb
                if grp == lb and r.next() % 2:
                    y = succ if (r.next() % 100000) / 100000.0 < min(1.0, ps * 0.6 + 0.3) else others[0]
                if r.next() % 25 == 0:
                    y = None
                g = grp if r.next() % 25 else None
                rows.append([y, g, y if r.next() % 3 else None])
            names = ["y", "g", "w"]
            level = LEVELS[int(r.next() % len(LEVELS))]
            p0t = P0S[int(r.next() % len(P0S))]
            if fi == 0:
                rows = [[succ, la, succ]]
            elif fi == 1:
                rows = [[succ, la, None]] * 5 + [[others[0], lb, None]] * 5
            elif fi == 2:
                rows = [[succ, la, succ]] * 7 + [[succ, lb, succ]] * 3
            elif fi == 3:
                p0t = "0.5"
                rows = [[succ, la, None]] * 5 + [[others[0], la, None]] * 5 + [[succ, lb, None]] * 4 + [[others[0], lb, None]] * 4
            n = len(rows)
            # --- one sample
            extra = ["--cols", "y,w", "--success", succ, "--p0", p0t, "--level", level]
            p = run(a.bin, tmp, names, rows, extra)
            if p.returncode != 0:
                keep = pathlib.Path(tempfile.gettempdir()) / ("gate_test4_failed_%d.csv" % fi)
                keep.write_text((tmp / "f.csv").read_text())
                print("gate_test4: exit %d on file %d: %s (the file is %s, args %r)" % (p.returncode, fi, p.stdout[:240], keep, extra))
                return 2
            doc = json.loads(p.stdout)["data"]
            want_level = float(Fraction(*[ot.level_parts(level)[0], 10 ** ot.level_parts(level)[1]]))
            if doc.get("test") != "prop" or doc.get("success") != succ or not gd.same(doc.get("p0"), float(Fraction(p0t))) or not gd.same(doc.get("level"), want_level) or doc.get("adjustment") != "none":
                bad.append((fi, "the header: %r" % {k: doc.get(k) for k in ("test", "success", "p0", "level", "adjustment")}))
            p0 = Fraction(p0t)
            for g, col in zip(doc["tests"], (0, 2)):
                cells = [row[col] for row in rows if row[col] is not None]
                xs = sum(1 for c in cells if c == succ)
                w = ot.prop_one(xs, len(cells), p0, level)
                check_one(bad, fi, "prop %s p0 %s level %s" % (names[col], p0t, level), g, w, n, p0, 10 ** 9, level, st)
                controls_one(st, g, w, p0, level)
            if n > 1:
                p2 = run(a.bin, tmp, names, shuffled(r, rows), extra)
                if p2.stdout != p.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
                p3 = run(a.bin, tmp, names, rows, ["--cols", "w,y", "--success", succ, "--p0", p0t, "--level", level])
                if p3.returncode != 0 or json.loads(p3.stdout)["data"]["tests"] != list(reversed(doc["tests"])):
                    bad.append((fi, "--cols in another order does not give the same entries in that order"))
            # --- two samples
            extra = ["--cols", "y,w", "--success", succ, "--group", "g", "--level", level]
            p = run(a.bin, tmp, names, rows, extra)
            gl = sorted(set(row[1] for row in rows if row[1] is not None))
            used_all = [row for row in rows if row[1] is not None]
            if len(gl) != 2:
                if p.returncode == 0:
                    bad.append((fi, "a group column of %d levels was not refused" % len(gl)))
                continue
            if p.returncode != 0:
                print("gate_test4: exit %d on file %d (two samples): %s" % (p.returncode, fi, p.stdout[:240]))
                return 2
            doc = json.loads(p.stdout)["data"]
            if doc.get("group") != "g" or "p0" in doc:
                bad.append((fi, "the header of two samples: %r" % {k: doc.get(k) for k in ("group", "p0")}))
            for g, col in zip(doc["tests"], (0, 2)):
                cnt = {}
                for row in used_all:
                    if row[col] is not None:
                        c = cnt.setdefault(row[1], [0, 0])
                        c[0] += 1
                        c[1] += 1 if row[col] == succ else 0
                (n1, x1), (n2, x2) = cnt.get(gl[0], [0, 0]), cnt.get(gl[1], [0, 0])
                w = ot.prop_two(x1, n1, x2, n2, level)
                check_two(bad, fi, "prop2 %s level %s" % (names[col], level), g, w, n, n1, x1, n2, x2, gl[0], gl[1], st)
                controls_two(st, g, w, n1, x1, n2, x2)
            if n > 1:
                p2 = run(a.bin, tmp, names, shuffled(r, rows), extra)
                if p2.stdout != p.stdout:
                    bad.append((fi, "two samples: the rows shuffled give different bytes"))
                swap = {la: "zz", lb: "aa"}
                srows = [[row[0], swap.get(row[1], row[1]), row[2]] for row in rows]
                p4 = run(a.bin, tmp, names, srows, extra)
                if p4.returncode == 0:
                    for g0, g4 in zip(doc["tests"], json.loads(p4.stdout)["data"]["tests"]):
                        ok = g0.get("p_value") == g4.get("p_value") and g0.get("reason") == g4.get("reason")
                        if g0["statistic"] is not None:
                            ok = ok and gd.same(-g0["statistic"] if g0["statistic"] != 0 else 0.0, g4["statistic"])
                        if g0.get("effect_size") is not None:
                            v0, v4 = g0["effect_size"]["value"], g4["effect_size"]["value"]
                            ok = ok and abs(v0 + v4) <= 1e-15 * max(1.0, abs(v0))
                        iv0, iv4 = g0.get("interval"), g4.get("interval")
                        if iv0 is not None:
                            ok = ok and gd.same(-iv0["estimate"] if iv0["estimate"] != 0 else 0.0, iv4["estimate"]) and abs(iv0["low"] + iv4["high"]) <= 1e-15 and abs(iv0["high"] + iv4["low"]) <= 1e-15
                        if not ok:
                            bad.append((fi, "the groups swapped do not negate the answer: %r against %r" % (g0, g4)))
    for k, (outside, total) in sorted(st.ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_test4 prop: %d files, %d entries compared, %d differ; p-values against mpmath: largest relative error %.3g over %d, at most %.3g of the bound; interval endpoints: largest absolute error %.3g (bound 1e-11); Cohen's h: largest absolute error %.3g (bound 1e-12); controls outside the bound: %s" % (
        a.files, st.checked, len(bad), st.p["max"], st.p["n"], st.p["ratio"], st.iv, st.h, ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(st.ctl.items()))))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
