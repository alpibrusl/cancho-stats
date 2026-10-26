#!/usr/bin/env python3
"""S4.2 gate R-D (calibration) and R-E (power, printed) of docs/rank.md section 5 for `stats test kruskal`:

    python3 scripts/gate_kruskal_calibration.py [--bin build/stats] [--replicates 4000] [--seed 1] [--workers 8]

Under a true null (k = 3 and 5 groups of one distribution, seeded: normal, and a lattice of 5 values and a binary one for heavy ties), for equal and unequal sizes, the rate at which `p < 0.05` lies inside the two-sided 99.9% binomial
interval for the number of replicates and the p-values are uniform (Kolmogorov-Smirnov, p > 0.001), for configurations with 5 or more values in every group (the chi-square approximation of H is printed, not gated, below:
it is known to be poor there). Where the p-values are on a lattice themselves (the lattice values, and normal values with a group below 10: H takes few values, the first run of (5,5,5) gave KS p 0.000 with 172 rejections inside the interval) only the rejection rate is judged. CONTROLS from the same replicates must be refused by the same test at some configuration, so
the gate is able to fail: chi-square on k degrees of freedom, H without the tie correction (on the lattices; the binary one is where it matters most, a tenth of the variation of a 5-value lattice: 2,000 replicates on the 5-value lattice did not refuse it, CI run of the first push), the one-sided p-value. Power against a planted shift (printed, not gated: there is no independent closed form for the
power of a rank test).

Exit 0: every test passes and every control is refused. Exit 1: a test fails. Exit 2: the harness failed.
"""
import argparse
import concurrent.futures
import json
import pathlib
import subprocess
import sys
import tempfile

import scipy.stats as ss

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
from gate_test_calibration import ALPHA_TEST, interval, normal  # noqa: E402

# (sizes, kind): normal values, or a lattice of 5 values
CONFIGS = [
    ((5, 5, 5), "normal"),
    ((10, 20, 40), "normal"),
    ((40, 20, 10), "normal"),
    ((5, 15, 40), "normal"),
    ((10, 10, 10, 10, 10), "normal"),
    ((3, 3, 3), "normal"),
    ((20, 20, 20), "lattice"),
    ((10, 30, 50), "lattice"),
    ((8, 8, 8, 8, 8), "lattice"),
    ((20, 20, 20), "binary"),
    ((12, 24, 36), "binary"),
]


def write(d, groups):
    rows = ["x,g"]
    for gi, vals in enumerate(groups):
        rows += ["%r,g%d" % (v, gi) for v in vals]
    (d / "f.csv").write_text("\n".join(rows) + "\n")


def one(binary, d, groups):
    write(d, groups)
    p = subprocess.run([binary, "test", "kruskal", "--root", str(d), "--cols", "x", "--group", "g", "f.csv"], capture_output=True, text=True, timeout=120)
    if p.returncode != 0:
        raise RuntimeError("exit %d: %s" % (p.returncode, p.stdout[:200]))
    return json.loads(p.stdout)["data"]["tests"][0]


def run_many(binary, jobs, workers):
    dirs = [pathlib.Path(tempfile.mkdtemp()) for _ in range(workers)]
    try:
        def work(k):
            return [one(binary, dirs[k], jobs[i]) for i in range(k, len(jobs), workers)]
        with concurrent.futures.ThreadPoolExecutor(workers) as ex:
            parts = list(ex.map(work, range(workers)))
        out = [None] * len(jobs)
        for k, part in enumerate(parts):
            for j, v in enumerate(part):
                out[k + j * workers] = v
        return out
    finally:
        for d in dirs:
            for f in d.glob("*"):
                f.unlink()
            d.rmdir()


def judge(ps, lo, hi, uniform=True):
    rej = sum(1 for p in ps if p < 0.05)
    ks = ss.kstest(ps, "uniform").pvalue if uniform else 1.0
    return lo <= rej <= hi and ks > ALPHA_TEST, rej, ks


def tie_term(groups):
    vals = [v for g in groups for v in g]
    n = len(vals)
    t3 = 0
    for v in set(vals):
        t = vals.count(v)
        t3 += t ** 3 - t
    return 1.0 - t3 / (n ** 3 - n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--replicates", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    reps = a.replicates
    lo, hi = interval(reps, 0.05)
    bad = []
    lines = []
    refused = {"chi-square-on-k": False, "no-tie-correction": False, "one-sided": False}
    try:
        for ci, (sizes, kind) in enumerate(CONFIGS):
            r = SM64(a.seed * 1000 + ci)
            if kind == "normal":
                jobs = [[[normal(r) for _ in range(n)] for n in sizes] for _ in range(reps)]
            else:
                top = 5 if kind == "lattice" else 2
                jobs = [[[float(int(r.next() % top)) for _ in range(n)] for n in sizes] for _ in range(reps)]
            res = run_many(a.bin, jobs, a.workers)
            k = len(sizes)
            ok_all = [x for x in res if x["p_value"] is not None]
            ps = [x["p_value"] for x in ok_all]
            gated = min(sizes) >= 5
            ok, rej, ks = judge(ps, lo, hi, uniform=(kind == "normal" and min(sizes) >= 10))
            lines.append("%s n=%s: rejected %d of %d (interval %d..%d)%s%s" % (kind, sizes, rej, len(ps), lo, hi, (", KS p %.3f" % ks) if kind == "normal" and min(sizes) >= 10 else "", ("" if ok else "  <-- FAILS") if gated else ("" if ok else "  (printed, not gated: a group below 5)")))
            if not ok and gated:
                bad.append("%s n=%s null calibration" % (kind, sizes))
            if gated:
                hs = [x["statistic"] for x in ok_all]
                okc, _, _ = judge([float(ss.chi2.sf(h, k)) for h in hs], lo, hi, uniform=(kind == "normal"))
                if not okc:
                    refused["chi-square-on-k"] = True
                okn, _, _ = judge([p / 2.0 for p in ps], lo, hi, uniform=(kind == "normal"))
                if not okn:
                    refused["one-sided"] = True
                if kind != "normal":
                    h0 = [h * tie_term(g) for h, g in zip(hs, [jobs[i] for i, x in enumerate(res) if x["p_value"] is not None])]
                    oku, _, _ = judge([float(ss.chi2.sf(h, k - 1)) for h in h0], lo, hi, uniform=False)
                    if not oku:
                        refused["no-tie-correction"] = True
        # power: printed
        sizes, shifts = (15, 20, 25), (0.0, 0.5, 1.0)
        pr = SM64(a.seed * 7919)
        jobs = [[[m + normal(pr) for _ in range(n)] for n, m in zip(sizes, shifts)] for _ in range(reps)]
        res = run_many(a.bin, jobs, a.workers)
        lines.append("planted shifts (0, 0.5, 1) sd, n=(15,20,25): rejected %d of %d (printed, not gated)" % (sum(1 for x in res if x["p_value"] < 0.05), reps))
    except Exception as e:
        print("gate_kruskal_calibration: harness failure: %s" % e)
        return 2
    for name, was in refused.items():
        lines.append("control %s: %s" % (name, "refused by the calibration test" if was else "NOT refused"))
        if not was:
            bad.append("control %s was not refused: the gate cannot tell it from the program" % name)
    print("gate_kruskal_calibration (%d replicates per configuration, seed %d):" % (reps, a.seed))
    for ln in lines:
        print("  " + ln)
    print("gate_kruskal_calibration: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
