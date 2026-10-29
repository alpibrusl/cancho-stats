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
import mpmath  # noqa: E402
import oracle_test as ot  # noqa: E402
from gate_test import numeric, parse  # noqa: E402
import gate_test2 as g2  # noqa: E402

STATS = ("mean", "median", "variance", "sd", "skewness", "kurtosis", "quantile", "trimmed-mean")
PROBS = ["0", "0.05", "0.25", "0.5", "0.9", "1", "0.333333333", "0.999"]
TRIMS = ["0", "0.1", "0.2", "0.25", "0.49", "0.333333333"]


def stat_args(stat, qt=6, prob="0.9", trim="0.2"):
    """the flags that go with a statistic"""
    if stat == "quantile":
        return ["--prob", prob, "--quantile-type", str(qt)]
    if stat == "trimmed-mean":
        return ["--trim", trim]
    return []
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


def is_num(x):
    return isinstance(x, float) and math.isfinite(x)


def safe_stat(name, vals, **kw):
    """the oracle's statistic: a float, inf where it is beyond the double, or the string "undefined" where it does not exist"""
    try:
        v = ot.boot_stat(name, vals, **kw)
    except (OverflowError, ValueError):
        return float("inf")
    except ZeroDivisionError:
        return "undefined"
    return float("inf") if v is None else v


def safe_reps(vals, name, reps, seed, kind=1, swap=False, **kw):
    n = len(vals)
    out = []
    for b in range(reps):
        picked = []
        for j in range(n):
            if swap:
                picked.append(vals[ot.boot_index(seed, kind, j, b, n)])
            else:
                picked.append(vals[ot.boot_index(seed, kind, b, j, n)])
        out.append(safe_stat(name, picked, **kw))
    return out


mpmath.mp.dps = 50


def mp_qnorm(p):
    return mpmath.sqrt(2) * mpmath.erfinv(2 * p - 1)


def bca_oracle(vs, stat, kw, est, reps_v, level, reps, zero_bias=False, zero_accel=False, ties_below=False):
    """the BCa quantities of docs/resample.md section 14 from the exact inputs, in mpmath at 50 digits: "undefined" where the section says bca-undefined, "extreme" where the program's double sums leave their range (judged
    only for not crashing), else a dict of bias, accel, low, high (the adjusted levels)"""
    n = len(vs)
    c = sum(1 for x in reps_v if x < est)
    e = sum(1 for x in reps_v if x == est)
    if ties_below:
        c, e = c + e, 0
    twice = 2 * c + e
    if twice <= 0 or twice >= 2 * reps:
        return "undefined"
    t = []
    for i in range(n):
        v = safe_stat(stat, vs[:i] + vs[i + 1:], **kw)
        if not is_num(v):
            return "undefined"
        t.append(v)
    ft = [Fraction(x) for x in t]
    m = sum(ft, Fraction(0)) / n
    s2 = sum(((m - x) ** 2 for x in ft), Fraction(0))
    s3 = sum(((m - x) ** 3 for x in ft), Fraction(0))
    if s2 == 0:
        return "undefined"
    if s2 < Fraction(1, 10 ** 280) or max(abs(m - x) for x in ft) > Fraction(10 ** 100):
        return "extreme"
    z0 = mp_qnorm(mpmath.mpf(twice) / (2 * reps))
    if zero_bias:
        z0 = mpmath.mpf(0)
    a = mpmath.mpf(s3.numerator) / s3.denominator / (6 * mpmath.power(mpmath.mpf(s2.numerator) / s2.denominator, mpmath.mpf(3) / 2))
    if zero_accel:
        a = mpmath.mpf(0)
    num, k = ot.level_parts(level)
    zl = mp_qnorm(mpmath.mpf(10 ** k - num) / (2 * 10 ** k))
    out = []
    for z in (zl, -zl):
        w = z0 + z
        d = 1 - a * w
        if d <= 0:
            return "undefined"
        out.append(mpmath.ncdf(z0 + w / d))
    return {"bias": z0, "accel": a, "low": out[0], "high": out[1]}


def close(got, want):
    return abs(got - float(want)) <= 1e-9 * abs(float(want)) + 1e-12


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
    und_seen = 0
    bca_undef = bca_extreme = bca_checked = 0

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
            stat = STATS[fi % len(STATS)]
            qt = 1 + int(r.next() % 9)
            prob = PROBS[int(r.next() % len(PROBS))]
            trim = TRIMS[int(r.next() % len(TRIMS))]
            kw = {}
            if stat == "quantile":
                kw = {"qt": qt, "p": Fraction(prob)}
            elif stat == "trimmed-mean":
                kw = {"trim": Fraction(trim)}
            seed = SEEDS[int(r.next() % len(SEEDS))]
            reps = REPS[int(r.next() % len(REPS))]
            if len(rows) > 60:
                reps = 100
            cols = ",".join(names[j] for j in num)
            method = "bca" if r.next() % 2 else "percentile"
            argv = ["ci", "boot", "--cols", cols, "--statistic", stat, "--seed", str(seed), "--replicates", str(reps), "--level", level] + stat_args(stat, qt, prob, trim) + (["--method", method] if method == "bca" or r.next() % 2 else [])
            p = run(a.bin, tmp, names, rows, argv)
            if p.returncode != 0:
                print("gate_ci_boot: exit %d on file %d: %s" % (p.returncode, fi, p.stdout[:300]))
                return 2
            d = json.loads(p.stdout)["data"]
            want_level = ot.rounded(Fraction(ot.level_parts(level)[0], 10 ** ot.level_parts(level)[1]))
            if (d.get("ci"), d.get("statistic"), d.get("method"), d.get("replicates"), d.get("seed"), d.get("rows")) != ("boot", stat, method, reps, seed, len(rows)) or d.get("level") != want_level:
                bad.append((fi, "the header: %r" % {k: d.get(k) for k in ("ci", "statistic", "method", "level", "replicates", "seed", "rows")}))
            if stat == "quantile" and (d.get("p"), d.get("type")) != (float(Fraction(prob)), qt):
                bad.append((fi, "the quantile header: %r" % {k: d.get(k) for k in ("p", "type")}))
            if stat == "trimmed-mean" and d.get("trim") != float(Fraction(trim)):
                bad.append((fi, "the trim header: %r" % d.get("trim")))
            if stat not in ("quantile", "trimmed-mean") and any(k in d for k in ("p", "type", "trim")):
                bad.append((fi, "a key of another statistic in the header: %r" % sorted(d)))
            lo_pos, hi_pos = ot.boot_positions(reps, level)
            if len(d["intervals"]) != len(num):
                bad.append((fi, "%d entries for %d columns" % (len(d["intervals"]), len(num))))
                continue
            for g, j in zip(d["intervals"], num):
                entries += 1
                vs = present(rows, j)
                if g["column"] != names[j] or g["n"] != {"n": len(vs), "missing": len(rows) - len(vs)}:
                    bad.append((fi, "%s: column/n %r" % (names[j], (g["column"], g["n"]))))
                est = safe_stat(stat, vs, **kw) if len(vs) >= (2 if stat in ("variance", "sd") else 1) else None
                if is_num(est):
                    if g["estimate"] != est:
                        bad.append((fi, "%s %s: estimate %r, want %r" % (stat, names[j], g["estimate"], est)))
                elif g["estimate"] is not None:
                    bad.append((fi, "%s %s: estimate %r where the oracle has none (%r)" % (stat, names[j], g["estimate"], est)))
                if len(vs) < 2:
                    nulls += 1
                    if g["interval"] is not None or g.get("reason") != "too-few-rows":
                        bad.append((fi, "%s: %d values: want a null interval with too-few-rows, got %r" % (names[j], len(vs), g)))
                    continue
                if est == "undefined":
                    nulls += 1
                    if g["interval"] is not None or g.get("reason") != "variance-zero":
                        bad.append((fi, "%s %s: constant column: want null with variance-zero, got %r" % (stat, names[j], g)))
                    continue
                reps_v = safe_reps(vs, stat, reps, seed, **kw)
                inf = not is_num(est) or any(isinstance(x, float) and not math.isfinite(x) for x in reps_v)
                und = sum(1 for x in reps_v if x == "undefined")
                if inf:
                    nulls += 1
                    if g["interval"] is not None or g.get("reason") != "statistic-beyond-double" or "undefined" in g:
                        bad.append((fi, "%s %s: want null with statistic-beyond-double, got %r" % (stat, names[j], g)))
                    continue
                if und:
                    nulls += 1
                    und_seen += 1
                    if g["interval"] is not None or g.get("reason") != "replicate-undefined" or g.get("undefined") != und:
                        bad.append((fi, "%s %s: want null with replicate-undefined and %d undefined, got %r" % (stat, names[j], und, g)))
                    continue
                srt = sorted(reps_v)
                iv = g["interval"]
                if method == "bca":
                    ora = bca_oracle(vs, stat, kw, est, reps_v, level, reps)
                    if ora == "undefined":
                        nulls += 1
                        bca_undef += 1
                        if iv is not None or g.get("reason") != "bca-undefined":
                            bad.append((fi, "%s %s: want null with bca-undefined, got %r" % (stat, names[j], g)))
                        continue
                    if ora == "extreme":
                        bca_extreme += 1
                        if iv is None and g.get("reason") != "bca-undefined":
                            bad.append((fi, "%s %s: extreme sums, null with %r" % (stat, names[j], g.get("reason"))))
                        continue
                    if iv is None:
                        bad.append((fi, "%s %s: a null interval where the BCa oracle has one (%s)" % (stat, names[j], g.get("reason"))))
                        continue
                    bca_checked += 1
                    adj = iv.get("adjusted", {})
                    if not (close(iv.get("bias", 1e9), ora["bias"]) and close(iv.get("acceleration", 1e9), ora["accel"]) and close(adj.get("low", 1e9), ora["low"]) and close(adj.get("high", 1e9), ora["high"])):
                        bad.append((fi, "%s %s: bias/acceleration/levels %r, want %s %s %s %s" % (stat, names[j], (iv.get("bias"), iv.get("acceleration"), adj), mpmath.nstr(ora["bias"], 17), mpmath.nstr(ora["accel"], 17), mpmath.nstr(ora["low"], 17), mpmath.nstr(ora["high"], 17))))
                        continue
                    pl = min(max(math.floor(reps * adj["low"]), 0), reps - 1)
                    ph = min(max(math.ceil(reps * adj["high"]) - 1, 0), reps - 1)
                    want = {"of": stat, "method": "bootstrap-bca", "level": want_level, "estimate": est, "low": srt[pl], "high": srt[ph], "positions": {"low": pl, "high": ph}, "bias": iv["bias"], "acceleration": iv["acceleration"], "adjusted": adj, "rounding": "double"}
                    if iv != want or list(iv) != list(want):
                        bad.append((fi, "%s %s seed %d R %d level %s: got %r want %r" % (stat, names[j], seed, reps, level, iv, want)))
                    for cname, kws in (("bca-no-bias", {"zero_bias": True}), ("bca-no-acceleration", {"zero_accel": True}), ("bca-ties-below", {"ties_below": True})):
                        if cname == "bca-ties-below" and not any(x == est for x in reps_v):
                            continue
                        alt = bca_oracle(vs, stat, kw, est, reps_v, level, reps, **kws)
                        if isinstance(alt, dict):
                            apl = min(max(math.floor(reps * float(alt["low"])), 0), reps - 1)
                            aph = min(max(math.ceil(reps * float(alt["high"])) - 1, 0), reps - 1)
                            control(cname, srt[apl] != srt[pl] or srt[aph] != srt[ph])
                    control("bca-round-positions", srt[min(max(round(reps * adj["low"]), 0), reps - 1)] != srt[pl] or srt[min(max(round(reps * adj["high"]) - 1, 0), reps - 1)] != srt[ph])
                    continue
                if iv is None:
                    bad.append((fi, "%s %s: a null interval with %d values (%s)" % (stat, names[j], len(vs), g.get("reason"))))
                    continue
                want = {"of": stat, "method": "bootstrap-percentile", "level": want_level, "estimate": est, "low": srt[lo_pos], "high": srt[hi_pos], "positions": {"low": lo_pos, "high": hi_pos}, "rounding": "double"}
                if iv != want or list(iv) != list(want):
                    bad.append((fi, "%s %s seed %d R %d level %s: got %r want %r" % (stat, names[j], seed, reps, level, iv, want)))
                # controls: other draws, other positions, other divisor
                alts = {
                    "another-seed": safe_reps(vs, stat, reps, seed ^ 1, **kw),
                    "kind-2": safe_reps(vs, stat, reps, seed, kind=2, **kw),
                    "indices-exchanged": safe_reps(vs, stat, reps, seed, swap=True, **kw),
                }
                if stat == "quantile" and qt != 7:
                    alts["quantile-type-7"] = safe_reps(vs, stat, reps, seed, **dict(kw, qt=7))
                if stat == "kurtosis":
                    alts["kurtosis-not-excess"] = [x + 3.0 for x in reps_v]
                if stat == "trimmed-mean":
                    n_ = len(vs)
                    alt_cut = int(Fraction(n_) * kw["trim"] + Fraction(1, 2))
                    if alt_cut != math.floor(n_ * kw["trim"]) and n_ - 2 * alt_cut >= 1:
                        rt = []
                        for b in range(reps):
                            pk_ = sorted(Fraction(vs[ot.boot_index(seed, 1, b, jj, n_)]) for jj in range(n_))
                            mid = pk_[alt_cut:n_ - alt_cut]
                            rt.append(ot.rounded(sum(mid, Fraction(0)) / len(mid)))
                        alts["trim-rounded"] = rt
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
            p1 = run(a.bin, tmp, ["x", "y"], [["3.5", "1"], ["", "2"], ["", "3"]], ["ci", "boot", "--cols", "x", "--statistic", st, "--seed", "5", "--replicates", "100"] + stat_args(st))
            g1 = json.loads(p1.stdout)["data"]["intervals"][0]
            want_est = 3.5 if st in ("mean", "median", "quantile", "trimmed-mean") else None
            if g1["estimate"] != want_est or g1["interval"] is not None or g1.get("reason") != "too-few-rows" or g1["n"] != {"n": 1, "missing": 2}:
                bad.append((-3, "one value, %s: %r" % (st, g1)))
        # small columns whose resamples are often constant: the rule for undefined replicates, whatever the random files reached
        for k, vs_ in enumerate(([1.0, 2.0, 5.0], [1.0, 1.0, 2.0], [0.5, 2.5, 4.0, 4.5])):
            for st in ("skewness", "kurtosis"):
                p1 = run(a.bin, tmp, ["x"], [[repr(v)] for v in vs_], ["ci", "boot", "--cols", "x", "--statistic", st, "--seed", str(3 + k), "--replicates", "100"])
                g1 = json.loads(p1.stdout)["data"]["intervals"][0]
                rv = safe_reps(vs_, st, 100, 3 + k)
                und1 = sum(1 for x in rv if x == "undefined")
                if und1:
                    ok1 = g1["interval"] is None and g1.get("reason") == "replicate-undefined" and g1.get("undefined") == und1
                else:
                    srt1 = sorted(rv)
                    lo1, hi1 = ot.boot_positions(100, "0.95")
                    ok1 = g1["interval"] is not None and (g1["interval"]["low"], g1["interval"]["high"]) == (srt1[lo1], srt1[hi1]) and "undefined" not in g1
                if not ok1:
                    bad.append((-4, "small column %r, %s: %r (the oracle has %d undefined)" % (vs_, st, g1, und1)))
                und_seen += 1 if und1 else 0
    if a.files >= 150 and und_seen == 0:
        bad.append((-1, "no column had an undefined replicate: the rule is not exercised"))
    for k, (differs, total) in sorted(ctl.items()):
        if total and differs * 10 < total:
            bad.append((-1, "control %s left the answer on only %d of %d cases: the gate cannot tell it from the program" % (k, differs, total)))
    print("gate_ci_boot: %d entries (%d null, %d with undefined replicates; BCa: %d checked, %d undefined as the oracle, %d extreme), %d problems; controls that differ from the program: %s" % (
        entries, nulls, und_seen, bca_checked, bca_undef, bca_extreme, len(bad), ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(ctl.items()))))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
