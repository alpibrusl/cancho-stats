#!/usr/bin/env python3
"""S5.1b gates I-A for `stats ci mean` and `stats ci diff` (docs/resample.md sections 3 and 9):

    python3 scripts/gate_ci_mean.py [--bin build/stats] [--files 200] [--seed 1]

The files of the gates of the tests (scripts/gate_test.py's shapes, gate_test2.py's groups of two labels with sizes down to 0, 1 and 2, missing cells, huge, tiny and constant values, levels from 0.5 to 0.999 and one of nine decimals)
run with `ci mean --cols` (against `test t1` on the same file), `ci mean --pairs` (against `test paired`) and `ci diff --cols --group` with `--variance welch` and `pooled` (against `test t2`). The `interval` object of every entry must be BYTE
EQUAL to the `interval` of the test, the counts equal, and the estimate equal to the interval's. Independently of the test: the endpoints of `ci mean` against scripts/oracle_test.py's `t1` (mpmath at 60 digits, the bound of T-C: absolute
1e-12 max(|estimate|, half-width), never below two units of 2^-1074) and of `ci diff` against `t2`. A null interval carries a `reason` (`too-few-rows` for fewer than 2 values in a sample). Also: the rows shuffled give the same bytes; the
columns named in another order give the same entries in that order; `--correct` is refused for `ci` (args.bad-value). CONTROLS from the same data must leave the bound on at least a tenth of the cases: the normal quantile in place of
Student's t, and the interval of the other variance rule (pooled for Welch and the reverse).

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
import oracle_test as ot  # noqa: E402
from gate_test import numeric, parse  # noqa: E402
import gate_test2 as g2  # noqa: E402
from gate_anova import shuffled  # noqa: E402


def run(binary, tmp, names, rows, argv):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary] + list(argv) + ["--root", str(tmp), "f.csv"], capture_output=True, text=True, timeout=300)


def vals_of(rows, j, gcol=None, label=None):
    out = []
    for rw in rows:
        if gcol is not None and rw[gcol] != label:
            continue
        v = parse(rw[j])
        if v is not None:
            out.append(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    ctl = {}
    counts = {"mean": 0, "paired": 0, "diff": 0}
    nulls = 0

    def control(name, outside):
        c = ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if outside else 0

    def same_entries(fi, tag, ci_doc, test_doc):
        ci_e = ci_doc["intervals"]
        te = test_doc["tests"]
        if len(ci_e) != len(te):
            bad.append((fi, "%s: %d entries against %d" % (tag, len(ci_e), len(te))))
            return
        for k, (g, h) in enumerate(zip(ci_e, te)):
            if g["interval"] != h["interval"]:
                bad.append((fi, "%s entry %d: ci %r against test %r" % (tag, k, g["interval"], h["interval"])))
            if g["interval"] is None:
                if g.get("reason") is None:
                    bad.append((fi, "%s entry %d: a null interval with no reason" % (tag, k)))
            elif g["estimate"] != g["interval"]["estimate"]:
                bad.append((fi, "%s entry %d: estimate %r against the interval's %r" % (tag, k, g["estimate"], g["interval"]["estimate"])))
            for key in ("n",):
                if key == "n" and "variance" not in g:
                    if g.get(key) != h.get(key):
                        bad.append((fi, "%s entry %d: n %r against %r" % (tag, k, g.get(key), h.get(key))))

    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            names, rows, gcol, labels, shapes, vidx, variance, level = g2.gen_t2(r, fi)
            num = [j for j in vidx if numeric([rw[j] for rw in rows])]
            if not num:
                continue
            cols = ",".join(names[j] for j in num)
            # ci mean against test t1
            pc = run(a.bin, tmp, names, rows, ["ci", "mean", "--cols", cols, "--level", level])
            pt = run(a.bin, tmp, names, rows, ["test", "t1", "--cols", cols, "--level", level])
            if pc.returncode != 0 or pt.returncode != 0:
                print("gate_ci_mean: exit %d / %d on file %d: %s | %s" % (pc.returncode, pt.returncode, fi, pc.stdout[:200], pt.stdout[:200]))
                return 2
            dc = json.loads(pc.stdout)["data"]
            dt = json.loads(pt.stdout)["data"]
            if dc.get("ci") != "mean" or not gd.same(dc.get("level"), ot.rounded(Fraction(*[ot.level_parts(level)[0], 10 ** ot.level_parts(level)[1]]))) or dc.get("rows") != len(rows):
                bad.append((fi, "the header of ci mean: %r" % {k: dc.get(k) for k in ("ci", "level", "rows")}))
            if "mu0" in dc or "adjustment" in dc or "tests" in dc:
                bad.append((fi, "the header of ci mean has a key of the test: %r" % sorted(dc)))
            same_entries(fi, "mean", dc, dt)
            counts["mean"] += len(dc["intervals"])
            for g, j in zip(dc["intervals"], num):
                vs = vals_of(rows, j)
                w = ot.t1(vs, 0.0, level)
                if g["interval"] is None:
                    nulls += 1
                    if len(vs) >= 2:
                        bad.append((fi, "mean %s: null interval with %d values" % (names[j], len(vs))))
                    continue
                iv = g["interval"]
                ok_lo = ot.endpoint_ok(iv["low"], w["low"], w["estimate"], w["half"])
                ok_hi = ot.endpoint_ok(iv["high"], w["high"], w["estimate"], w["half"])
                if not ok_lo or not ok_hi:
                    bad.append((fi, "mean %s level %s: want %s %s got %r %r" % (names[j], level, mpmath.nstr(w["low"], 17), mpmath.nstr(w["high"], 17), iv["low"], iv["high"])))
                if w["half"] is not None and w["half"] != 0 and len(vs) >= 3:
                    n = len(vs)
                    z = ot.z_quantile(level) if hasattr(ot, "z_quantile") else None
                    if z is not None:
                        half_z = mpmath.mpf(z) * (w["half"] / w["q"])
                        control("normal-for-t", not ot.endpoint_ok(float(w["estimate"] - float(half_z)), w["low"], w["estimate"], w["half"]))
            # the rows shuffled; the columns in another order
            if len(rows) > 1:
                p2 = run(a.bin, tmp, names, shuffled(r, rows), ["ci", "mean", "--cols", cols, "--level", level])
                if p2.stdout != pc.stdout:
                    bad.append((fi, "the rows shuffled give different bytes"))
                if len(num) > 1:
                    p3 = run(a.bin, tmp, names, rows, ["ci", "mean", "--cols", ",".join(reversed([names[j] for j in num])), "--level", level])
                    d3 = {x["column"]: x for x in json.loads(p3.stdout)["data"]["intervals"]}
                    if any(d3[x["column"]] != x for x in dc["intervals"]):
                        bad.append((fi, "--cols in another order changes an entry"))
            # ci diff against test t2, both variance rules
            gl = sorted(set(rw[gcol] for rw in rows if rw[gcol] != ""))
            for vr in ("welch", "pooled"):
                pd_ = run(a.bin, tmp, names, rows, ["ci", "diff", "--cols", cols, "--group", names[gcol], "--variance", vr, "--level", level])
                pt2 = run(a.bin, tmp, names, rows, ["test", "t2", "--cols", cols, "--group", names[gcol], "--variance", vr, "--level", level])
                if len(gl) != 2:
                    if pd_.returncode != pt2.returncode or json.loads(pd_.stdout).get("error", {}).get("rule") != json.loads(pt2.stdout).get("error", {}).get("rule"):
                        bad.append((fi, "diff with %d levels: ci exit %d, test exit %d" % (len(gl), pd_.returncode, pt2.returncode)))
                    continue
                if pd_.returncode != 0 or pt2.returncode != 0:
                    bad.append((fi, "diff %s: exit %d / %d %s" % (vr, pd_.returncode, pt2.returncode, pd_.stdout[:120])))
                    continue
                dd = json.loads(pd_.stdout)["data"]
                dt2 = json.loads(pt2.stdout)["data"]
                if dd.get("ci") != "diff" or dd.get("group") != names[gcol] or dd.get("variance") != vr or dd.get("levels") != dt2.get("levels"):
                    bad.append((fi, "the header of ci diff: %r" % {k: dd.get(k) for k in ("ci", "group", "variance", "levels")}))
                if "adjustment" in dd or "tests" in dd:
                    bad.append((fi, "the header of ci diff has a key of the test: %r" % sorted(dd)))
                same_entries(fi, "diff " + vr, dd, dt2)
                counts["diff"] += len(dd["intervals"])
                for g, j in zip(dd["intervals"], num):
                    va = vals_of(rows, j, gcol, gl[0])
                    vb = vals_of(rows, j, gcol, gl[1])
                    w = ot.t2(va, vb, vr == "pooled", level)
                    if g["interval"] is None:
                        nulls += 1
                        if len(va) >= 2 and len(vb) >= 2 and w["low"] is not None:
                            bad.append((fi, "diff %s: null interval where the oracle has one" % names[j]))
                        continue
                    if w["low"] is None:
                        continue
                    iv = g["interval"]
                    q = w.get("q", 0) or 0
                    if not ot.endpoint_ok(iv["low"], w["low"], w["estimate"], w["half"], q) or not ot.endpoint_ok(iv["high"], w["high"], w["estimate"], w["half"], q):
                        bad.append((fi, "diff %s %s: want %s %s got %r %r" % (vr, names[j], mpmath.nstr(w["low"], 17), mpmath.nstr(w["high"], 17), iv["low"], iv["high"])))
                    other = ot.t2(va, vb, vr != "pooled", level)
                    if other["low"] is not None and w["half"] != 0:
                        control("other-variance-rule", not ot.endpoint_ok(float(other["low"]), w["low"], w["estimate"], w["half"], q))
            # ci mean --pairs against test paired
        for fi in range(a.files):
            names, rows, pairs, shapes, level = g2.gen_paired(r, fi)
            if not pairs:
                continue
            ptxt = ",".join("%s:%s" % pr for pr in pairs)
            pc = run(a.bin, tmp, names, rows, ["ci", "mean", "--pairs", ptxt, "--level", level])
            pt = run(a.bin, tmp, names, rows, ["test", "paired", "--pairs", ptxt, "--level", level])
            if pc.returncode != pt.returncode:
                bad.append((1000 + fi, "paired: ci exit %d, test exit %d: %s" % (pc.returncode, pt.returncode, pc.stdout[:100])))
                continue
            if pc.returncode != 0:
                continue
            dc = json.loads(pc.stdout)["data"]
            dt = json.loads(pt.stdout)["data"]
            if dc.get("ci") != "mean" or "mu0" in dc or "adjustment" in dc:
                bad.append((1000 + fi, "the header of ci mean --pairs: %r" % sorted(dc)))
            same_entries(1000 + fi, "paired", dc, dt)
            counts["paired"] += len(dc["intervals"])
        # columns of one value, of one value among empty cells, and the key of a pair
        for k, (cells, what) in enumerate((([["3.5"]], "one value"), ([["3.5"], [""], [""]], "one value among empty cells"), ([["2"], ["2"], ["2"]], "every value equal"), ([["1"], ["2"]], "two values"))):
            rows1 = [[c[0], c[0]] for c in cells]
            nm = ["x", "y"]
            pc = run(a.bin, tmp, nm, rows1, ["ci", "mean", "--cols", "x"])
            pt = run(a.bin, tmp, nm, rows1, ["test", "t1", "--cols", "x"])
            if pc.returncode != 0 or pt.returncode != 0:
                bad.append((2000 + k, "%s: exit %d / %d" % (what, pc.returncode, pt.returncode)))
                continue
            g = json.loads(pc.stdout)["data"]["intervals"][0]
            h = json.loads(pt.stdout)["data"]["tests"][0]
            vs = [float(c[0]) for c in cells if c[0] != ""]
            want_est = sum(vs) / len(vs)
            if g["interval"] != h["interval"] or g.get("estimate") != want_est or g["n"] != h["n"] or g["column"] != "x":
                bad.append((2000 + k, "%s: ci %r against test %r" % (what, g, h)))
            if len(vs) < 2 and (g["interval"] is not None or g.get("reason") != "too-few-rows"):
                bad.append((2000 + k, "%s: want a null interval with too-few-rows, got %r" % (what, g)))
        rows2 = [["1", "2"], ["2", "4"], ["4", "5"]]
        pc = run(a.bin, tmp, ["a", "b"], rows2, ["ci", "mean", "--pairs", "a:b"])
        g = json.loads(pc.stdout)["data"]["intervals"][0]
        if g.get("pair") != ["a", "b"] or "column" in g:
            bad.append((2100, "the pair entry has the key(s) %r" % sorted(g)))
        # --correct is not for ci
        names, rows, gcol, labels, shapes, vidx, variance, level = g2.gen_t2(SM64(3), 40)
        p = run(a.bin, tmp, names, rows, ["ci", "mean", "--cols", names[vidx[0]], "--correct", "holm"])
        doc = json.loads(p.stdout) if p.stdout else {}
        if p.returncode != 2 or doc.get("error", {}).get("rule") != "args.bad-value":
            bad.append((-2, "--correct with ci: want exit 2 args.bad-value, got %d %r" % (p.returncode, p.stdout[:120])))
    for k, (outside, total) in sorted(ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s left the bound on only %d of %d cases: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_ci_mean: entries compared with the tests: mean %d, paired %d, diff %d (%d null intervals), %d problems; controls outside the bound: %s" % (
        counts["mean"], counts["paired"], counts["diff"], nulls, len(bad), ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(ctl.items()))))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
