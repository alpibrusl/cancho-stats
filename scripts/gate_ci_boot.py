#!/usr/bin/env python3
"""S5.2a gates I-C and I-D for `stats ci boot` (docs/resample.md sections 4 and 12):

    python3 scripts/gate_ci_boot.py [--bin build/stats] [--files 60] [--seed 1]

The files of the gates of the tests (scripts/gate_test2.py's shapes: missing cells, huge, tiny and constant values, columns of 2 to 200 rows, levels from 0.5 to 0.999 and one of nine decimals) run with `ci boot` for the four
statistics, a seed below 2^62 (the largest among them) and 100 to 400 replicates. An independent implementation (scripts/oracle_test.py: SplitMix64 from the specification, checked against its published outputs; the same
indices; the statistic of each resample as a Fraction rounded once) recomputes ALL the replicates and the program's `low` and `high` must be bit-equal to the replicates at the printed positions of the sorted list, the positions
must follow the rule in exact integers, the estimate must be the statistic of the data, and a column that cannot have an interval is null with its reason. Also: the same bytes on a second run and with the columns named in
another order (the entries in that order). CONTROLS from the same data must leave the answer on at least a tenth of the cases: another seed, the kind 2 in place of 1, the replicate index and the draw index exchanged, the
positions floor(R (1 + level) / 2) for the high one, and the divisor n - 1 for a mean.

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

STATS = ("mean", "median", "variance", "sd")
SEEDS = [0, 1, 7, 20240229, 4611686018427387903, 1234567890123456789, 2 ** 40 + 5]
REPS = [100, 120, 160, 200, 240, 400, 800]


def run(binary, tmp, names, rows, argv):
    text = ",".join(names) + "\n" + "\n".join(",".join(r) for r in rows) + "\n"
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary] + list(argv) + ["--root", str(tmp), "f.csv"], capture_output=True, text=True, timeout=600)


def present(rows, j):
    out = []
    for rw in rows:
        v = parse(rw[j])
        if v is not None:
            out.append(v)
    return out


def finite_all(xs):
    return all(isinstance(x, float) and math.isfinite(x) for x in xs)


def safe_stat(name, vals):
    try:
        v = ot.boot_stat(name, vals)
    except (OverflowError, ValueError):
        return float("inf")
    return float("inf") if v is None else v


def safe_reps(vals, name, reps, seed, kind=1, swap=False):
    n = len(vals)
    out = []
    for b in range(reps):
        picked = []
        for j in range(n):
            if swap:
                picked.append(vals[ot.boot_index(seed, kind, j, b, n)])
            else:
                picked.append(vals[ot.boot_index(seed, kind, b, j, n)])
        out.append(safe_stat(name, picked))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=60)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    # the SplitMix64 of the oracle against the published outputs (seed 0, the first three)
    if [ot.sm_mix((k + 1) * ot.GOLD & ot.MASK64) for k in range(3)] != [0xE220A8397B1DCDAF, 0x6E789E6AA1B965F4, 0x06C45D188009454F]:
        print("gate_ci_boot: the oracle's SplitMix64 is not the published one")
        return 2
    r = SM64(a.seed)
    bad = []
    ctl = {}
    entries = 0
    nulls = 0

    def control(name, differs):
        c = ctl.setdefault(name, [0, 0])
        c[1] += 1
        c[0] += 1 if differs else 0

    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            names, rows, gcol, labels, shapes, vidx, variance, level = g2.gen_t2(r, fi)
            num = [j for j in vidx if numeric([rw[j] for rw in rows])]
            if not num:
                continue
            stat = STATS[fi % 4]
            seed = SEEDS[int(r.next() % len(SEEDS))]
            reps = REPS[int(r.next() % len(REPS))]
            if len(rows) > 60:
                reps = 100
            cols = ",".join(names[j] for j in num)
            argv = ["ci", "boot", "--cols", cols, "--statistic", stat, "--seed", str(seed), "--replicates", str(reps), "--level", level]
            p = run(a.bin, tmp, names, rows, argv)
            if p.returncode != 0:
                print("gate_ci_boot: exit %d on file %d: %s" % (p.returncode, fi, p.stdout[:300]))
                return 2
            d = json.loads(p.stdout)["data"]
            want_level = ot.rounded(Fraction(ot.level_parts(level)[0], 10 ** ot.level_parts(level)[1]))
            if (d.get("ci"), d.get("statistic"), d.get("method"), d.get("replicates"), d.get("seed"), d.get("rows")) != ("boot", stat, "percentile", reps, seed, len(rows)) or d.get("level") != want_level:
                bad.append((fi, "the header: %r" % {k: d.get(k) for k in ("ci", "statistic", "method", "level", "replicates", "seed", "rows")}))
            lo_pos, hi_pos = ot.boot_positions(reps, level)
            if len(d["intervals"]) != len(num):
                bad.append((fi, "%d entries for %d columns" % (len(d["intervals"]), len(num))))
                continue
            for g, j in zip(d["intervals"], num):
                entries += 1
                vs = present(rows, j)
                if g["column"] != names[j] or g["n"] != {"n": len(vs), "missing": len(rows) - len(vs)}:
                    bad.append((fi, "%s: column/n %r" % (names[j], (g["column"], g["n"]))))
                est = safe_stat(stat, vs) if len(vs) >= (2 if stat in ("variance", "sd") else 1) else None
                if est is None and g["estimate"] is not None:
                    bad.append((fi, "%s %s: estimate %r of %d values" % (stat, names[j], g["estimate"], len(vs))))
                if est is not None and math.isfinite(est):
                    if g["estimate"] != est:
                        bad.append((fi, "%s %s: estimate %r, want %r" % (stat, names[j], g["estimate"], est)))
                elif est is not None and g["estimate"] is not None:
                    bad.append((fi, "%s %s: estimate %r where the oracle's is beyond the double" % (stat, names[j], g["estimate"])))
                if len(vs) < 2:
                    nulls += 1
                    if g["interval"] is not None or g.get("reason") != "too-few-rows":
                        bad.append((fi, "%s: %d values: want a null interval with too-few-rows, got %r" % (names[j], len(vs), g)))
                    continue
                reps_v = safe_reps(vs, stat, reps, seed)
                if not finite_all(reps_v) or est is None or not math.isfinite(est):
                    nulls += 1
                    if g["interval"] is not None or g.get("reason") != "statistic-beyond-double":
                        bad.append((fi, "%s %s: want null with statistic-beyond-double, got %r" % (stat, names[j], g)))
                    continue
                srt = sorted(reps_v)
                iv = g["interval"]
                if iv is None:
                    bad.append((fi, "%s %s: a null interval with %d values (%s)" % (stat, names[j], len(vs), g.get("reason"))))
                    continue
                want = {"of": stat, "method": "bootstrap-percentile", "level": want_level, "estimate": est, "low": srt[lo_pos], "high": srt[hi_pos], "positions": {"low": lo_pos, "high": hi_pos}, "rounding": "double"}
                if iv != want or list(iv) != list(want):
                    bad.append((fi, "%s %s seed %d R %d level %s: got %r want %r" % (stat, names[j], seed, reps, level, iv, want)))
                # controls: other draws, other positions, other divisor
                alts = {
                    "another-seed": safe_reps(vs, stat, reps, seed ^ 1),
                    "kind-2": safe_reps(vs, stat, reps, seed, kind=2),
                    "indices-exchanged": safe_reps(vs, stat, reps, seed, swap=True),
                }
                for k, rv in alts.items():
                    if finite_all(rv):
                        s2 = sorted(rv)
                        control(k, s2[lo_pos] != srt[lo_pos] or s2[hi_pos] != srt[hi_pos])
                num_l, kd = ot.level_parts(level)
                alt_hi = min(reps - 1, reps * (10 ** kd + num_l) // (2 * 10 ** kd))
                control("floor-high-position", srt[alt_hi] != srt[hi_pos])
                if stat == "mean":
                    n = len(vs)
                    rm = []
                    for b in range(reps):
                        s = sum((Fraction(vs[ot.boot_index(seed, 1, b, jj, n)]) for jj in range(n)), Fraction(0))
                        rm.append(ot.rounded(s / (n - 1)))
                    s3 = sorted(rm)
                    control("n-minus-1", s3[lo_pos] != srt[lo_pos] or s3[hi_pos] != srt[hi_pos])
            # same bytes again; the columns in another order
            p2 = run(a.bin, tmp, names, rows, argv)
            if p2.stdout != p.stdout:
                bad.append((fi, "a second run gives other bytes"))
            if len(num) > 1:
                rev = list(reversed([names[j] for j in num]))
                argv3 = list(argv)
                argv3[3] = ",".join(rev)
                p3 = run(a.bin, tmp, names, rows, argv3)
                d3 = {x["column"]: x for x in json.loads(p3.stdout)["data"]["intervals"]}
                if [x["column"] for x in json.loads(p3.stdout)["data"]["intervals"]] != [names[j] for j in num] or any(d3[x["column"]] != x for x in d["intervals"]):
                    bad.append((fi, "--cols in another order changes an entry or the order of the entries"))
        # a column of one value: the estimate of the data for the mean and the median, null for the variance and the sd; never an interval
        for st in STATS:
            p1 = run(a.bin, tmp, ["x", "y"], [["3.5", "1"], ["", "2"], ["", "3"]], ["ci", "boot", "--cols", "x", "--statistic", st, "--seed", "5", "--replicates", "100"])
            g1 = json.loads(p1.stdout)["data"]["intervals"][0]
            want_est = 3.5 if st in ("mean", "median") else None
            if g1["estimate"] != want_est or g1["interval"] is not None or g1.get("reason") != "too-few-rows" or g1["n"] != {"n": 1, "missing": 2}:
                bad.append((-3, "one value, %s: %r" % (st, g1)))
    for k, (differs, total) in sorted(ctl.items()):
        if total and differs * 10 < total:
            bad.append((-1, "control %s left the answer on only %d of %d cases: the gate cannot tell it from the program" % (k, differs, total)))
    print("gate_ci_boot: %d entries (%d null), %d problems; controls that differ from the program: %s" % (
        entries, nulls, len(bad), ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(ctl.items()))))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
