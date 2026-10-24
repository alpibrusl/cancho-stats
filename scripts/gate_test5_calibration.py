#!/usr/bin/env python3
"""S3.5 gate T-G (docs/test.md section 12), the error rates of the corrections:

    python3 scripts/gate_test5_calibration.py [--bin build/stats] [--replicates 4000] [--seed 1] [--workers 8]

Twenty independent one-sample t tests (n = 30, normal data; `stats test t1 --cols c0,..,c19 --correct ..`), seeded. The p-values of independent true nulls are uniform and independent, so the rates have exact values:
  * global null (all 20 true): the probability of at least one adjusted p-value below 0.05 is 1 - (1 - 0.05 / 20)^20 = 0.04889 for Bonferroni and for Holm (the smallest p-value decides) and exactly 0.05 for
    Benjamini-Hochberg (its false discovery rate under independence is pi0 alpha, and with every null true a discovery is a false one);
  * half of the nulls true (columns with a shift of 10 standard errors of the mean in the other ten, so that Holm has rejected them all before it reaches a true null): Holm's family-wise error is 1 - (1 - 0.05 / 10)^10 = 0.04897 and Bonferroni's 1 - (1 - 0.05 / 20)^10 = 0.02471 and the false discovery rate of BH is
    pi0 alpha = 0.025, which the mean false discovery proportion of the replicates must show within three standard errors.
The count of replicates with a false rejection must lie inside the two-sided 99.9% binomial interval around the exact probability. CONTROLS: the p-values uncorrected (their family-wise rate is 1 - 0.95^20 = 0.642)
must be refused by the same test, and so must Holm's rate read as BH's (0.04889 against 0.05 is not distinguishable: not used); the uncorrected control is the one that has to fail.

Exit 0: every test passes and the control is refused. Exit 1: a test fails. Exit 2: the harness failed.
"""
import argparse
import concurrent.futures
import json
import math
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
from gate_test_calibration import interval, normal  # noqa: E402

M = 20
N = 30
SHIFT = 10.0 / math.sqrt(N)   # ten standard errors of the mean, in units of the standard deviation: the ten false nulls are rejected before the true ones are reached


def one(binary, d, job):
    cols = job
    (d / "f.csv").write_text(",".join("c%d" % j for j in range(M)) + "\n" + "\n".join(",".join(repr(cols[j][i]) for j in range(M)) for i in range(N)) + "\n")
    out = {}
    for m in ("bonferroni", "holm", "bh"):
        p = subprocess.run([binary, "test", "t1", "--cols", ",".join("c%d" % j for j in range(M)), "--mu0", "0", "--correct", m, "--root", str(d), "f.csv"], capture_output=True, text=True, timeout=120)
        if p.returncode != 0:
            raise RuntimeError("exit %d: %s" % (p.returncode, p.stdout[:300]))
        doc = json.loads(p.stdout)["data"]
        out[m] = doc["family"]["p_adjusted"]
        out["p"] = [t["p_value"] for t in doc["tests"]]
    return out


def run_many(binary, jobs, workers):
    dirs = [pathlib.Path(tempfile.mkdtemp()) for _ in range(workers)]
    try:
        def work(k):
            return [one(binary, dirs[k], jobs[i]) for i in range(k, len(jobs), workers)]
        with concurrent.futures.ThreadPoolExecutor(workers) as ex:
            parts = list(ex.map(work, range(workers)))
        res = [None] * len(jobs)
        for k, part in enumerate(parts):
            for j, v in enumerate(part):
                res[k + j * workers] = v
        return res
    finally:
        for d in dirs:
            for f in d.glob("*"):
                f.unlink()
            d.rmdir()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--replicates", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    reps = a.replicates
    bad, lines = [], []
    refused = False
    try:
        for name, true_nulls in (("global null", M), ("half of the nulls true", M // 2)):
            r = SM64(a.seed * 10 + true_nulls)
            jobs = []
            for _ in range(reps):
                cols = []
                for j in range(M):
                    shift = 0.0 if j < true_nulls else SHIFT
                    cols.append([shift + normal(r) for _ in range(N)])
                jobs.append(cols)
            res = run_many(a.bin, jobs, a.workers)
            def any_false(vals):
                return sum(1 for v in vals if any(x < 0.05 for x in v[:true_nulls]))
            fw = 1 - (1 - 0.05 / true_nulls) ** true_nulls
            for m in ("bonferroni", "holm"):
                if m == "bonferroni":
                    fw_m = 1 - (1 - 0.05 / M) ** true_nulls
                else:
                    fw_m = fw
                cnt = any_false([x[m] for x in res])
                lo, hi = interval(reps, fw_m)
                ok = lo <= cnt <= hi
                lines.append("%s, %s: family-wise error %d of %d (exact %.5f, interval %d..%d)%s" % (name, m, cnt, reps, fw_m, lo, hi, "" if ok else "  <-- FAILS"))
                if not ok:
                    bad.append("%s %s family-wise error" % (name, m))
            # BH: false discovery proportion
            fdp = []
            for x in res:
                adj = x["bh"]
                rej = [i for i, v in enumerate(adj) if v < 0.05]
                fdp.append(sum(1 for i in rej if i < true_nulls) / len(rej) if rej else 0.0)
            mean = sum(fdp) / reps
            sd = math.sqrt(sum((v - mean) ** 2 for v in fdp) / (reps - 1))
            want = (true_nulls / M) * 0.05
            se = sd / math.sqrt(reps)
            ok = abs(mean - want) <= 3 * se
            lines.append("%s, bh: mean false discovery proportion %.5f (exact %.5f, three standard errors %.5f)%s" % (name, mean, want, 3 * se, "" if ok else "  <-- FAILS"))
            if not ok:
                bad.append("%s bh false discovery rate" % name)
            # the control: the uncorrected p-values
            cnt = any_false([x["p"] for x in res])
            lo, hi = interval(reps, fw)
            if not (lo <= cnt <= hi):
                refused = True
            lines.append("control no correction, %s: family-wise error %d of %d (inside the interval of the corrected rate: %s)" % (name, cnt, reps, "yes" if lo <= cnt <= hi else "no"))
    except Exception as e:
        print("gate_test5_calibration: harness failure: %s" % e)
        return 2
    lines.append("control no correction: %s" % ("refused by the test" if refused else "NOT refused"))
    if not refused:
        bad.append("the control without correction was not refused: the gate cannot tell it from the program")
    print("gate_test5_calibration (%d replicates, %d columns of %d values, seed %d):" % (reps, M, N, a.seed))
    for ln in lines:
        print("  " + ln)
    print("gate_test5_calibration: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
