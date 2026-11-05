#!/usr/bin/env python3
"""S6.3 gates G-C1 to G-C4 for `stats cycle detect` (docs/estimation.md sections 1.2, 3 and 7):

    python3 scripts/gate_cycledetect.py [--bin build/stats] [--cases 40] [--seed 1]

scripts/oracle_cycle.py is an independent implementation written from the document (numpy for the least squares, SplitMix64 for the shuffles of the document's rule): on seeded multi-series files, with and without a cycle, the
decision (`found`, its reason), the p-value (a count, so EXACTLY), the grid (its size and ends), the gain of every period of the grid and the best one, and for a found cycle the refined period, the amplitude and the gain must equal
the oracle's: the gains and the amplitude to a relative 1e-9, the refined period to a relative 1e-6 (a ternary search on a flat maximum cannot resolve a period better than the square root of the rounding). The same bytes with the rows
in another order. Series shorter than 3 (order + 3) days are dropped and counted. CONTROLS that must leave the answer on at least a tenth of the cases: another stream kind for the shuffles, one shuffle fewer, a level of one more
degree. Every refusal is reached with its tag, and "nothing found" is an answer (exit 0).

Exit 0: every number is the expected one. Exit 1: a difference. Exit 2: the harness failed.
"""
import argparse
import json
import math
import pathlib
import subprocess
import sys
import tempfile

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
import oracle_cycle as oc  # noqa: E402


def run(binary, tmp, text, argv):
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "cycle", "detect"] + list(argv) + ["--root", str(tmp), "f.csv"], capture_output=True, text=True, timeout=900)


def parse(r):
    d = json.loads(r.stdout)
    return d["ok"], d.get("data"), d.get("error")


def rel(a, b, floor=1e-12):
    return abs(a - b) / (abs(b) + floor)


def make(rng):
    nser = rng.pick([1, 2, 4, 8])
    cyc = rng.pick([0, 1, 1, 1])
    period = rng.pick([30.0, 40.0, 55.0, 90.0, 120.0])
    amp = rng.pick([0.04, 0.07, 0.1, 0.15, 0.5]) if cyc else 0.0
    n = rng.pick([120, 200, 270, 330])
    use_w = rng.next() % 3 == 0
    gaps = rng.next() % 4 == 0
    series = {}
    for k in range(nser):
        days, d = [], 0
        for i in range(n if k % 5 != 4 else 10):
            days.append(d)
            d += 1 + (1 if gaps and rng.next() % 5 == 0 else 0)
        ph = rng.u() * 6.283
        ys = [round(amp * math.sin(2 * math.pi * x / period + ph) + 0.002 * x + (rng.u() - 0.5) * 1.2, 6) for x in days]
        ws = [round(0.3 + rng.u() * 2, 4) if use_w else 1.0 for _ in days]
        series["ser%02d" % k] = (days, ys, ws)
    return series, use_w


def text_of(series, use_w, reverse=False):
    rows = [(k, d, y, w) for k, (ds, ys, ws) in series.items() for d, y, w in zip(ds, ys, ws)]
    if reverse:
        rows.reverse()
    return "k,day,y" + (",w" if use_w else "") + "\n" + "\n".join(",".join(str(v) for v in (r if use_w else r[:3])) for r in rows) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--cases", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-controls", action="store_true", help="skip the controls (the mutants run uses this: a control is a statement about the gate, not about the build)")
    a = ap.parse_args()
    rng = SM64(a.seed)
    bad = []
    worst = {"gain": 0.0, "amp": 0.0, "period": 0.0, "grid": 0.0}
    differs = {"kind": 0, "fewer": 0, "order": 0}
    found_n = 0
    done = 0
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        for k in range(a.cases):
            series, use_w = make(rng)
            order = rng.pick([2, 2, 2, 1, 3])
            shuffles = rng.pick([99, 199, 299])
            seed = rng.pick([0, 1, 7, 20240229, 2 ** 40 + 5])
            minp = rng.pick([21.0, 21.0, 14.0])
            alpha = rng.pick([0.004, 0.05, 0.2, 0.5])
            if k % 3 == 1:
                # a series of exactly 3 (order + 3) days is used, one day fewer is dropped
                need_ = 3 * (order + 3)
                series["ser97"] = (list(range(need_)), [round(math.sin(i / 2.0) + 0.1 * (i % 3), 6) for i in range(need_)], [1.0] * need_)
                series["ser98"] = (list(range(need_ - 1)), [round(math.cos(i / 2.0), 6) for i in range(need_ - 1)], [1.0] * (need_ - 1))
                if use_w:
                    series["ser97"] = (series["ser97"][0], series["ser97"][1], [round(0.5 + 0.1 * (i % 4), 4) for i in range(need_)])
            argv = ["--day", "day", "--y", "y", "--series", "k", "--seed", str(seed), "--order", str(order), "--shuffles", str(shuffles), "--min-period", repr(minp), "--alpha", repr(alpha), "--grid", "yes"] + (["--weight", "w"] if use_w else [])
            r = run(a.bin, tmp, text_of(series, use_w), argv)
            ok, data, err = parse(r)
            if not ok:
                bad.append(("case %d refused" % k, err))
                continue
            r2 = run(a.bin, tmp, text_of(series, use_w, reverse=True), argv)
            if r2.stdout != r.stdout:
                bad.append(("case %d order" % k,))
            keys = sorted(series)
            need = 3 * (order + 3)
            used = [(np.array(series[q][0], float), np.array(series[q][1], float), np.array(series[q][2], float)) for q in keys if sum(1 for w in series[q][2] if w > 0) >= need]
            want = oc.detect(used, seed, order=order, min_period=minp, permutations=shuffles, alpha=alpha)
            if data["series_used"] != len(used) or data["series_dropped"] != len(keys) - len(used):
                bad.append(("case %d series counts" % k, data["series_used"], data["series_dropped"], len(used)))
            if data["found"] != want["found"]:
                bad.append(("case %d found" % k, data["found"], want["found"], data.get("reason"), want.get("reason")))
                continue
            if not data["found"]:
                if data["reason"] != want["reason"]:
                    bad.append(("case %d reason" % k, data["reason"], want["reason"]))
            if want.get("reason") != "too-short":
                if data["p"] != want["p"]:
                    bad.append(("case %d p" % k, data["p"], want["p"]))
                gains = data["gains"]
                if len(gains) != len(want["gains"]) or data["grid"]["periods"] != len(want["grid"]):
                    bad.append(("case %d grid size" % k, len(gains), len(want["gains"])))
                else:
                    scale = max(abs(g) for g in want["gains"]) + 1e-12
                    gg = max(abs(x - y) for x, y in zip(gains, want["gains"])) / scale
                    worst["gain"] = max(worst["gain"], gg)
                    if gg > 1e-9:
                        bad.append(("case %d gains" % k, gg))
                    if rel(data["grid"]["shortest"], want["grid"][-1]) > 1e-12 or rel(data["grid"]["longest"], want["grid"][0]) > 1e-12:
                        bad.append(("case %d grid ends" % k, data["grid"]))
                    if rel(data["best_on_grid"]["period"], want["best_period"]) > 1e-12 or rel(data["best_on_grid"]["gain"], want["best_gain"]) > 1e-9:
                        bad.append(("case %d best" % k, data["best_on_grid"], want["best_period"], want["best_gain"]))
            if data["found"]:
                found_n += 1
                dp, da, dg = rel(data["period"], want["period"]), rel(data["amplitude"], want["amplitude"]), rel(data["gain"], want["gain"])
                worst["period"], worst["amp"] = max(worst["period"], dp), max(worst["amp"], da)
                if dp > 1e-6 or da > 1e-9 or dg > 1e-9:
                    bad.append(("case %d refined" % k, data["period"], want["period"], data["amplitude"], want["amplitude"]))
            # controls (on the first cases only: the oracle is slow)
            if k < 14 and not a.no_controls and want.get("reason") != "too-short":
                for name, kw in (("kind", {"kind": 5}), ("fewer", {"permutations": shuffles - 1}), ("order", {"order": order + 1 if order < 3 else order - 1})):
                    args_ = dict(order=order, min_period=minp, permutations=shuffles, alpha=alpha)
                    args_.update(kw)
                    w2 = oc.detect(used, seed, **args_)
                    if w2.get("p") != want.get("p") or w2.get("best_gain") != want.get("best_gain") or w2["found"] != want["found"]:
                        differs[name] += 1
            done += 1
        for name, count in differs.items():
            if not a.no_controls and count * 10 < min(14, a.cases):
                bad.append(("control %s left the answer on %d of %d cases" % (name, count, min(14, a.cases)),))
        if found_n == 0 or found_n == done:
            bad.append(("the cases found %d cycles of %d: the gate cannot tell a decision from its opposite" % (found_n, done),))
        # ---- nothing found is an answer; too short is an answer
        rng2 = SM64(5)
        noise = {"s": (list(range(150)), [round((rng2.u() - 0.5), 6) for _ in range(150)], [1.0] * 150)}
        r = run(a.bin, tmp, text_of(noise, False), ["--day", "day", "--y", "y", "--series", "k", "--seed", "1"])
        ok, data, err = parse(r)
        if not ok or r.returncode != 0 or data["found"] or data["reason"] not in ("p-above-alpha", "edge"):
            bad.append(("noise", ok, data and data.get("reason"), r.returncode))
        short = {"s": (list(range(30)), [round(math.sin(i / 4), 6) for i in range(30)], [1.0] * 30)}
        r = run(a.bin, tmp, text_of(short, False), ["--day", "day", "--y", "y", "--series", "k", "--seed", "1"])
        ok, data, err = parse(r)
        if not ok or data["found"] or data["reason"] != "too-short" or r.returncode != 0:
            bad.append(("too short", ok, data and data.get("reason")))
        tiny = {"s": (list(range(10)), [1.0] * 10, [1.0] * 10)}
        r = run(a.bin, tmp, text_of(tiny, False), ["--day", "day", "--y", "y", "--series", "k", "--seed", "1"])
        ok, data, err = parse(r)
        if not ok or data["found"] or data["series_dropped"] != 1 or data["series_used"] != 0 or data["reason"] != "too-short":
            bad.append(("dropped series", data, err))
        # ---- refusals
        base = "k,day,y\n" + "\n".join("s,%d,%r" % (i, math.sin(i)) for i in range(60)) + "\n"
        argv0 = ["--day", "day", "--y", "y", "--series", "k", "--seed", "1"]
        refusals = [
            ("stats.seed-required", base, ["--day", "day", "--y", "y", "--series", "k"]),
            ("args.bad-value", base, ["--day", "day", "--y", "y", "--series", "k", "--seed", "x"]),
            ("stats.test-columns", base, ["--day", "day", "--y", "y", "--seed", "1"]),
            ("args.bad-value", base, argv0 + ["--min-period", "1"]),
            ("args.bad-value", base, argv0 + ["--min-cycles", "0.5"]),
            ("args.bad-value", base, argv0 + ["--alpha", "0"]),
            ("args.bad-value", base, argv0 + ["--alpha", "1"]),
            ("args.bad-value", base, argv0 + ["--edge", "1.5"]),
            ("args.bad-value", base, argv0 + ["--edge", "0"]),
            ("args.bad-value", base, argv0 + ["--shuffles", "98"]),
            ("args.bad-value", base, argv0 + ["--shuffles", "10000"]),
            ("args.bad-value", base, argv0 + ["--oversample", "0"]),
            ("args.bad-value", base, argv0 + ["--oversample", "17"]),
            ("args.bad-value", base, argv0 + ["--order", "4"]),
            ("stats.no-such-column", base, ["--day", "nope", "--y", "y", "--series", "k", "--seed", "1"]),
            ("stats.cycle.duplicate-day", base + "s,3,1\n", argv0),
            ("stats.cycle.day-not-integer", "k,day,y\ns,1.5,2\ns,2,3\n", argv0),
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
    print("cycle detect: %d cases (%d found a cycle); worst relative gap gains %.2e period %.2e amplitude %.2e; controls differed on %s of %d" % (done, found_n, worst["gain"], worst["period"], worst["amp"], differs, min(14, a.cases)))
    if bad:
        for b in bad[:20]:
            print("FAIL", *b)
        print("%d differences" % len(bad))
        return 1
    print("gate_cycledetect: ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
