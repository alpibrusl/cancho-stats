#!/usr/bin/env python3
"""G-C2 to G-C5 of docs/estimation.md section 7, measured (not asserted) for `stats cycle detect`:

    python3 scripts/measure_cycle.py [--bin build/stats] [--runs 10]

G-C2 power: a cycle of 40, 55, 90 or 120 days, three amplitudes, 270 days of noise uniform in +-0.6, one series, `runs` seeds each: the share found and the largest error of the period.
G-C3 seed sensitivity: 30 weak-cycle series, 20 seeds each: the share of series whose decision is not the same under all the seeds.
G-C4 false positives: `fp` series (2000 by default) of noise, with growth and with a step at day 135, one series of 270 days, default alpha: the share found.
G-C5 cost: 20 series of 270 days, the minimum of 5 runs.
Prints the numbers; the exit status is 0 (it is a measurement, the floors are in the document).
"""
import argparse
import json
import math
import pathlib
import subprocess
import sys
import tempfile
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402


def detect(binary, tmp, rows, seed, extra=()):
    (tmp / "f.csv").write_text("k,day,y\n" + "\n".join("%s,%d,%r" % r for r in rows) + "\n")
    r = subprocess.run([binary, "cycle", "detect", "--day", "day", "--y", "y", "--series", "k", "--seed", str(seed), "--root", str(tmp), "f.csv"] + list(extra), capture_output=True, text=True, timeout=600)
    return json.loads(r.stdout)["data"]


def series(rng, key, n=270, period=None, amp=0.0, trend=0.0, step=0.0):
    ph = rng.u() * 6.283
    out = []
    for i in range(n):
        y = (rng.u() - 0.5) * 1.2 + trend * i + (step if i >= n // 2 else 0.0)
        if period:
            y += amp * math.sin(2 * math.pi * i / period + ph)
        out.append((key, i, round(y, 6)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--runs", type=int, default=10)
    ap.add_argument("--fp", type=int, default=2000, help="series per false-positive setting")
    a = ap.parse_args()
    rng = SM64(11)
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        print("G-C2 power (270 days, noise +-0.6, one series, %d seeds each)" % a.runs)
        for period in (40, 55, 90, 120):
            for amp in (0.1, 0.2, 0.4):
                found, err = 0, 0.0
                for s in range(a.runs):
                    d = detect(a.bin, tmp, series(rng, "s", period=period, amp=amp), s + 1)
                    if d["found"]:
                        found += 1
                        err = max(err, abs(d["period"] - period) / period)
                print("  period %3d amplitude %.1f: found %2d of %d, largest period error %.1f%%" % (period, amp, found, a.runs, 100 * err))
        flips = 0
        for k in range(30):
            rows = series(rng, "s", period=90, amp=0.1)
            ds = {detect(a.bin, tmp, rows, s)["found"] for s in range(1, 21)}
            flips += len(ds) > 1
        print("G-C3 seed sensitivity: %d of 30 weak-cycle series change decision with the seed (20 seeds each)" % flips)
        for name, kw in (("noise", {}), ("growth", {"trend": 0.004}), ("step", {"step": 0.6})):
            f = 0
            for k in range(a.fp):
                f += detect(a.bin, tmp, series(rng, "s", **kw), k + 1)["found"]
            p = f / a.fp
            print("G-C4 false positives, %s: %d of %d found (%.2f%%, binomial standard error %.2f points; alpha 0.004)" % (name, f, a.fp, 100 * p, 100 * math.sqrt(max(p * (1 - p), 1e-9) / a.fp)))
        rows = []
        for k in range(20):
            rows += series(rng, "s%02d" % k, period=90, amp=0.3)
        ts = []
        for _ in range(5):
            t0 = time.perf_counter()
            detect(a.bin, tmp, rows, 7)
            ts.append(time.perf_counter() - t0)
        print("G-C5 cost: 20 series of 270 days, minimum of 5 runs %.2f s (median %.2f s)" % (min(ts), sorted(ts)[2]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
