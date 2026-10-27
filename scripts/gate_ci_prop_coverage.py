#!/usr/bin/env python3
"""S5.1 gate I-B (docs/resample.md section 6) for the intervals of `stats ci prop`, one sample:

    python3 scripts/gate_ci_prop_coverage.py [--bin build/stats] [--level 0.95]

The coverage of an interval method at (n, p) is a sum over the number of successes x of the binomial probabilities of the x whose interval contains p: no simulation, so there is no sampling noise to hide behind. The intervals of
every x = 0 .. n (n = 5, 10, 20, 40, 80) come from the PROGRAM (one run for each x, `--method wilson` and `--method clopper-pearson`), the probabilities from mpmath at 30 digits, on a grid of 600 values of p (a grid of irrational
steps, so that no p sits on an endpoint, plus p near 0 and 1). Properties, each fixed before the first run and corrected in place in docs/resample.md if a run shows it false:
  (1) Clopper-Pearson's coverage is at least the level, at every p of the grid and every n (the method's guarantee; a tolerance of 1e-9 for the rounding of the endpoints);
  (2) Wilson's mean absolute distance of the coverage from the level, over the grid, is smaller than the Wald interval's at every n (the reason Wilson is the default);
  (3) the Wald interval (computed here from x, n and the normal quantile; the CONTROL) has a coverage below the level minus 0.05 at some p for every n: the gate can tell a bad interval from a good one.
The coverages are printed (minimum, mean, and the p of the minimum) so that the table can go in the document.

Exit 0: the properties hold. Exit 1: one fails. Exit 2: the harness failed.
"""
import argparse
import json
import math
import pathlib
import subprocess
import sys
import tempfile

import mpmath

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import oracle_test as ot  # noqa: E402

mpmath.mp.dps = 30


def intervals(binary, tmp, n, method, level):
    out = []
    for x in range(n + 1):
        rows = ["y"] * x + ["n"] * (n - x)
        (tmp / "f.csv").write_text("o\n" + "\n".join(rows) + "\n")
        p = subprocess.run([binary, "ci", "prop", "--root", str(tmp), "--cols", "o", "--success", "y", "--level", level, "--method", method, "f.csv"], capture_output=True, text=True, timeout=60)
        if p.returncode != 0:
            raise RuntimeError("exit %d: %s" % (p.returncode, p.stdout[:200]))
        iv = json.loads(p.stdout)["data"]["intervals"][0]["interval"]
        out.append((iv["low"], iv["high"]))
    return out


def coverage(n, ivs, p):
    q = mpmath.mpf(p)
    tot = mpmath.mpf(0)
    for x, (lo, hi) in enumerate(ivs):
        if lo <= p <= hi:
            tot += mpmath.binomial(n, x) * q ** x * (1 - q) ** (n - x)
    return float(tot)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--level", default="0.95")
    a = ap.parse_args()
    level = float(a.level)
    num, k = ot.level_parts(a.level)
    z = float(mpmath.sqrt(2) * mpmath.erfinv(mpmath.mpf(num) / 10 ** k))
    grid = [0.0005, 0.001, 0.005, 0.995, 0.999, 0.9995] + [((i + 1) * 0.6180339887498949) % 1.0 * 0.998 + 0.001 for i in range(594)]
    bad = []
    lines = []
    try:
        with tempfile.TemporaryDirectory() as t:
            tmp = pathlib.Path(t)
            for n in (5, 10, 20, 40, 80):
                wil = intervals(a.bin, tmp, n, "wilson", a.level)
                cpi = intervals(a.bin, tmp, n, "clopper-pearson", a.level)
                wald = []
                for x in range(n + 1):
                    ph = x / n
                    h = z * math.sqrt(ph * (1 - ph) / n)
                    wald.append((ph - h, ph + h))
                cw = [coverage(n, wil, p) for p in grid]
                cc = [coverage(n, cpi, p) for p in grid]
                cd = [coverage(n, wald, p) for p in grid]
                mw = sum(abs(c - level) for c in cw) / len(grid)
                md = sum(abs(c - level) for c in cd) / len(grid)
                lines.append("n=%3d  CP min %.4f mean %.4f | Wilson min %.4f mean %.4f mean|dist| %.4f | Wald min %.4f (p %.4f) mean %.4f mean|dist| %.4f" % (
                    n, min(cc), sum(cc) / len(cc), min(cw), sum(cw) / len(cw), mw, min(cd), grid[cd.index(min(cd))], sum(cd) / len(cd), md))
                if min(cc) < level - 1e-9:
                    bad.append("(1) n=%d: Clopper-Pearson coverage %.6f at p=%.5f is below the level" % (n, min(cc), grid[cc.index(min(cc))]))
                if not mw < md:
                    bad.append("(2) n=%d: Wilson's mean distance %.4f is not below Wald's %.4f" % (n, mw, md))
                if not min(cd) < level - 0.05:
                    bad.append("(3) n=%d: the Wald control has minimum coverage %.4f, not below the level minus 0.05: the gate cannot tell" % (n, min(cd)))
    except Exception as e:
        print("gate_ci_prop_coverage: harness failure: %s" % e)
        return 2
    print("gate_ci_prop_coverage (level %s, 600 values of p):" % a.level)
    for ln in lines:
        print("  " + ln)
    print("gate_ci_prop_coverage: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
