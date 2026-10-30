#!/usr/bin/env python3
"""S5.4a gates for `stats resample perm` (docs/resample.md section 15):

    python3 scripts/gate_resample_perm.py [--bin build/stats] [--files 60] [--seed 1]

The files of the gates of the tests (scripts/gate_test2.py's two labels, missing cells, huge, tiny and constant values) run with `resample perm --cols --group` in the Monte Carlo mode and, where the arrangements are few, with
`--permutations all`. An independent implementation (scripts/oracle_test.py: the shuffle from the specification, SplitMix64 checked against its published outputs, the statistic and the comparison as Fractions) recomputes the count
of extreme replicates, which must equal the printed `extreme` EXACTLY; `p` must be the correctly rounded (1 + c) / (N + 1) (c / C for the exact test), `statistic` the mean difference rounded once, `mc_error` within relative 1e-12.
Also: the same bytes on a second run and with the columns in another order. CONTROLS from the same data must leave the count or the p on at least a tenth of the cases: another seed, kind 1 in place of 2, `>` in place of `>=`, a one-sided
count, c / N in place of (1 + c) / (N + 1), the first n_b values as group a*.

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
import oracle_test as ot  # noqa: E402
from gate_test import numeric, parse  # noqa: E402
import gate_test2 as g2  # noqa: E402

SEEDS = [0, 1, 7, 20240229, 4611686018427387903, 1234567890123456789, 2 ** 40 + 5]
PERMS = [100, 120, 150, 200, 250]


def run(binary, tmp, names, rows, argv):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "resample", "perm"] + list(argv) + ["--root", str(tmp), "f.csv"], capture_output=True, text=True, timeout=600)


def vals_of(rows, j, gcol, label):
    out = []
    for rw in rows:
        if rw[gcol] != label:
            continue
        v = parse(rw[j])
        if v is not None:
            out.append(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=60)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if [ot.sm_mix((k + 1) * ot.GOLD & ot.MASK64) for k in range(3)] != [0xE220A8397B1DCDAF, 0x6E789E6AA1B965F4, 0x06C45D188009454F]:
        print("gate_resample_perm: the oracle's SplitMix64 is not the published one")
        return 2
    r = SM64(a.seed)
    bad = []
    ctl = {}
    entries = exact_entries = nulls = 0

    def control(name, differs):
        c = ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if differs else 0

    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            names, rows, gcol, labels, shapes, vidx, variance, level = g2.gen_t2(r, fi)
            num = [j for j in vidx if numeric([rw[j] for rw in rows])]
            gl = sorted(set(rw[gcol] for rw in rows if rw[gcol] != ""))
            if not num or len(gl) != 2:
                continue
            seed = SEEDS[int(r.next() % len(SEEDS))]
            perms = PERMS[int(r.next() % len(PERMS))]
            if len(rows) > 60:
                perms = 100
            cols = ",".join(names[j] for j in num)
            for mode in ("mc", "all"):
                if mode == "all":
                    # the exact test only where every column has few arrangements
                    from math import comb
                    if any(comb(len(vals_of(rows, j, gcol, gl[0])) + len(vals_of(rows, j, gcol, gl[1])), len(vals_of(rows, j, gcol, gl[0]))) > 4000 for j in num):
                        continue
                    argv = ["--cols", cols, "--group", names[gcol], "--permutations", "all"]
                else:
                    argv = ["--cols", cols, "--group", names[gcol], "--seed", str(seed), "--permutations", str(perms)]
                p = run(a.bin, tmp, names, rows, argv)
                if p.returncode != 0:
                    print("gate_resample_perm: exit %d on file %d: %s" % (p.returncode, fi, p.stdout[:300]))
                    return 2
                d = json.loads(p.stdout)["data"]
                if (d.get("resample"), d.get("statistic"), d.get("group"), d.get("levels"), d.get("rows"), d.get("p_method")) != ("perm", "mean-diff", names[gcol], {"a": gl[0], "b": gl[1]}, len(rows), "exact" if mode == "all" else "monte-carlo"):
                    bad.append((fi, "the header: %r" % {k: d.get(k) for k in ("resample", "statistic", "group", "levels", "rows", "p_method")}))
                if mode == "mc" and d.get("seed") != seed or mode == "all" and "seed" in d:
                    bad.append((fi, "the seed key: %r" % d.get("seed")))
                if len(d["tests"]) != len(num):
                    bad.append((fi, "%d entries for %d columns" % (len(d["tests"]), len(num))))
                    continue
                for g, j in zip(d["tests"], num):
                    entries += 1
                    va, vb = vals_of(rows, j, gcol, gl[0]), vals_of(rows, j, gcol, gl[1])
                    inrows = sum(1 for rw in rows if rw[gcol] in gl)
                    want_n = {"a": len(va), "b": len(vb), "missing": inrows - len(va) - len(vb)}
                    if g["column"] != names[j] or g["n"] != want_n:
                        bad.append((fi, "%s: column/n %r want %r" % (names[j], (g["column"], g["n"]), want_n)))
                    if not va or not vb:
                        nulls += 1
                        if g["p"] is not None or g.get("reason") != "too-few-rows" or g["statistic"] is not None:
                            bad.append((fi, "%s: %d/%d values: want a null entry with too-few-rows, got %r" % (names[j], len(va), len(vb), g)))
                        continue
                    if mode == "mc":
                        c, tt = ot.perm_mc(va, vb, seed, perms)
                        total = perms
                        want_p = ot.rounded(Fraction(1 + c, total + 1))
                    else:
                        exact_entries += 1
                        c, total, tt = ot.perm_all(va, vb)
                        want_p = ot.rounded(Fraction(c, total))
                    want_t = ot.rounded(tt)
                    if want_t is None:
                        continue
                    if g["extreme"] != c or g["permutations"] != total or g["p"] != want_p or g["statistic"] != want_t:
                        bad.append((fi, "%s %s seed %d: got extreme %r permutations %r p %r T %r, want %d %d %r %r" % (mode, names[j], seed, g["extreme"], g["permutations"], g["p"], g["statistic"], c, total, want_p, want_t)))
                        continue
                    if mode == "mc":
                        mce = math.sqrt(want_p * (1 - want_p) / perms)
                        if abs(g.get("mc_error", 1e9) - mce) > 1e-12 * mce + 1e-300:
                            bad.append((fi, "%s: mc_error %r want %r" % (names[j], g.get("mc_error"), mce)))
                    elif "mc_error" in g:
                        bad.append((fi, "%s: mc_error in an exact test" % names[j]))
                    # controls
                    if mode == "mc":
                        for cname, kw in (("another-seed", {"seed": seed ^ 1}), ("kind-1", {"kind": 1}), ("strict", {"strict": True}), ("one-sided", {"one_sided": True}), ("first-nb", {"first_nb": True})):
                            args = dict(seed=seed, perms=perms)
                            args.update(kw)
                            c2, _ = ot.perm_mc(va, vb, **args)
                            control(cname, c2 != c)
                        control("c-over-N", ot.rounded(Fraction(c, perms)) != want_p)
                    else:
                        c2, _, _ = ot.perm_all(va, vb, strict=True)
                        control("exact-strict", c2 != c)
                        c3, _, _ = ot.perm_all(va, vb, one_sided=True)
                        control("exact-one-sided", c3 != c)
                p2 = run(a.bin, tmp, names, rows, argv)
                if p2.stdout != p.stdout:
                    bad.append((fi, "a second run gives other bytes"))
                if len(num) > 1:
                    argv3 = list(argv)
                    argv3[1] = ",".join(reversed([names[j] for j in num]))
                    p3 = run(a.bin, tmp, names, rows, argv3)
                    t3 = json.loads(p3.stdout)["data"]["tests"]
                    if [x["column"] for x in t3] != list(reversed([names[j] for j in num])) and [x["column"] for x in t3] != [names[j] for j in num]:
                        bad.append((fi, "--cols in another order changes the entries"))
                    by = {x["column"]: x for x in t3}
                    if any(by[x["column"]] != x for x in d["tests"]):
                        bad.append((fi, "--cols in another order changes an entry"))
        # fixed: two arrangements, every value equal, a negative zero
        for k, (va_, vb_) in enumerate((([1.0], [2.0]), ([3.0, 3.0], [3.0, 3.0, 3.0]), ([0.0, -0.0], [0.0]))):
            rows1 = [[repr(v), "a"] for v in va_] + [[repr(v), "b"] for v in vb_]
            p1 = run(a.bin, tmp, ["x", "g"], rows1, ["--cols", "x", "--group", "g", "--seed", "4", "--permutations", "100"])
            g1 = json.loads(p1.stdout)["data"]["tests"][0]
            c, tt = ot.perm_mc(va_, vb_, 4, 100)
            if g1["extreme"] != c or g1["p"] != ot.rounded(Fraction(1 + c, 101)):
                bad.append((-4 - k, "fixed case %r %r: %r want extreme %d" % (va_, vb_, g1, c)))
        # a group with no value in the column: a null entry, not a refusal
        p1 = run(a.bin, tmp, ["x", "g"], [["1", "a"], ["2", "a"], ["", "b"], ["", "b"]], ["--cols", "x", "--group", "g", "--seed", "4", "--permutations", "100"])
        g1 = json.loads(p1.stdout)["data"]["tests"][0]
        if g1["p"] is not None or g1.get("reason") != "too-few-rows" or g1["n"] != {"a": 2, "b": 0, "missing": 2}:
            bad.append((-8, "a group with no value: %r" % g1))
    for k, (differs, total) in sorted(ctl.items()):
        if total and differs * 10 < total:
            bad.append((-1, "control %s left the answer on only %d of %d cases: the gate cannot tell it from the program" % (k, differs, total)))
    print("gate_resample_perm: %d entries (%d exact, %d null), %d problems; controls that differ from the program: %s" % (
        entries, exact_entries, nulls, len(bad), ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(ctl.items()))))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
