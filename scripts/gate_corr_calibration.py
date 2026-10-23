#!/usr/bin/env python3
"""S2.2 gates C-E (calibration) and C-F (power) of docs/corr.md section 5, for the Pearson p-value of `stats corr`.

    python3 scripts/gate_corr_calibration.py [--bin build/stats] [--replicates 4000] [--power-replicates 4000] [--seed 1] [--workers 8]

C-E. Under a true null (two independent normal columns, seeded) and for n = 5, 10, 30, 100, 1000: the rate at which the program's p-value is below 0.05 lies inside the two-sided 99.9% binomial
interval for the number of replicates, and the p-values are uniform (Kolmogorov-Smirnov against U(0,1), p > 0.001). CONTROLS, computed from the same replicates and the same r, with the program's p-value
replaced by a wrong one, must be refused by the same two tests at some n, so the gate is able to fail: a one-sided tail (p / 2), the normal distribution in place of Student's t, and n in place of
n - 2 degrees of freedom.

C-F. With rho = 0.5 and n = 50 (bivariate normal, seeded) the rate at which p < 0.05 equals the exact power: the integral over |r| > r_crit of Fisher's density of the sample correlation (mpmath, 30
digits), r_crit from the 0.975 quantile of Student's t on 48 degrees of freedom; inside the 99.9% binomial interval. The density is first checked to integrate to 1.

A rate compared is a count of p-values below 0.05 printed as decimal text; nothing here judges a digit of r or p (that is scripts/gate_corr.py). Exit 0: every test passes and every control is refused.
Exit 1: a test fails or a control is not refused. Exit 2: the harness failed.
"""
import argparse
import concurrent.futures
import json
import math
import pathlib
import subprocess
import sys
import tempfile

import mpmath
import scipy.stats as ss

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402

mpmath.mp.dps = 30
NS = (5, 10, 30, 100, 1000)
ALPHA_TEST = 0.001


def normal(r):
    """a standard normal draw from the seeded stream (Box-Muller; the library log and cos only shape the test data, never a judged value)"""
    u1 = max(r.u(), 1e-300)
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * r.u())


def sample(r, n, rho):
    xs = [normal(r) for _ in range(n)]
    zs = [normal(r) for _ in range(n)]
    ys = [rho * x + math.sqrt(1.0 - rho * rho) * z for x, z in zip(xs, zs)]
    return xs, ys


METHOD = "pearson"


def corr_p(binary, d, xs, ys):
    f = d / "f.csv"
    f.write_text("x,y\n" + "\n".join("%r,%r" % (a, b) for a, b in zip(xs, ys)) + "\n")
    p = subprocess.run([binary, "corr", "--root", str(d), "--method", METHOD, "f.csv"], capture_output=True, text=True, timeout=120)
    if p.returncode != 0:
        raise RuntimeError("exit %d: %s" % (p.returncode, p.stdout[:200]))
    pair = json.loads(p.stdout)["data"]["pairs"][0]
    return pair[{"spearman": "rho", "kendall": "tau"}.get(METHOD, "r")], pair["p_value"]


def run_many(binary, jobs, workers):
    """jobs: list of (xs, ys); answers [(r, p)] in order"""
    dirs = [pathlib.Path(tempfile.mkdtemp()) for _ in range(workers)]
    try:
        def work(k):
            return [corr_p(binary, dirs[k], *jobs[i]) for i in range(k, len(jobs), workers)]
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


def binom_interval(reps, p=0.05):
    lo = int(ss.binom.ppf(ALPHA_TEST / 2, reps, p))
    hi = int(ss.binom.ppf(1 - ALPHA_TEST / 2, reps, p))
    return lo, hi


def fisher_density(r, rho, n):
    """the density of the sample correlation of a bivariate normal sample of size n (Fisher 1915; Johnson and Kotz): n >= 4"""
    r = mpmath.mpf(r)
    rho = mpmath.mpf(rho)
    c = (n - 2) * mpmath.gamma(n - 1) * (1 - rho ** 2) ** (mpmath.mpf(n - 1) / 2) / (mpmath.sqrt(2 * mpmath.pi) * mpmath.gamma(n - mpmath.mpf(1) / 2))
    return c * (1 - r ** 2) ** (mpmath.mpf(n - 4) / 2) / (1 - rho * r) ** (n - mpmath.mpf(3) / 2) * mpmath.hyp2f1(0.5, 0.5, n - mpmath.mpf(1) / 2, (1 + rho * r) / 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--replicates", type=int, default=4000)
    ap.add_argument("--power-replicates", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--method", choices=["pearson", "spearman", "kendall"], default="pearson")
    a = ap.parse_args()
    global METHOD
    METHOD = a.method
    bad = []
    lo, hi = binom_interval(a.replicates)
    refused = {"one-sided": False, "normal": False, "df-n": False}
    if a.method == "kendall":
        refused = {"one-sided": False}  # the normal-for-t and n-degrees-of-freedom controls are controls of a t test, not of the normal approximation of tau
    lines = []
    try:
        for n in NS:
            r = SM64(a.seed * 1000 + n)
            jobs = [sample(r, n, 0.0) for _ in range(a.replicates)]
            res = run_many(a.bin, jobs, a.workers)
            ps = [p for _, p in res]
            rej = sum(1 for p in ps if p < 0.05)
            ks = ss.kstest(ps, "uniform").pvalue
            gated = a.method == "pearson" or n >= 30
            # the p-values of Kendall's tau at n = 30 are on a lattice of about 400 values and the normal approximation to them is visibly not uniform to a KS test at these replicates (KS p 0.001, 0.000, 0.023 for
            # seeds 1 to 3, the rejection rate at 0.05 inside its interval every time): the rate is gated from n = 30, the uniformity from n = 100
            ks_gated = ks > ALPHA_TEST or (a.method == "kendall" and n < 100)
            ok = lo <= rej <= hi and ks_gated
            lines.append("n=%d: rejected %d of %d (interval %d..%d), KS p %.3f%s" % (n, rej, a.replicates, lo, hi, ks, ("" if gated else "  (printed, not gated: the approximation of the p-value is rough at small n)") if ok else ("  <-- FAILS" if gated else "  (outside; printed, not gated)")))
            if not ok and gated:
                bad.append("n=%d null calibration" % n)
            # controls: the same r, a wrong p
            for name in refused:
                cp = []
                for (rr, p) in res:
                    if rr is None or abs(rr) >= 1.0:
                        cp.append(p)
                        continue
                    t = rr * math.sqrt((n - 2) / (1.0 - rr * rr))
                    if name == "one-sided":
                        cp.append(p / 2.0)
                    elif name == "normal":
                        cp.append(2.0 * ss.norm.sf(abs(t)))
                    else:
                        cp.append(2.0 * ss.t.sf(abs(t), n))
                crej = sum(1 for p in cp if p < 0.05)
                cks = ss.kstest(cp, "uniform").pvalue
                if not (lo <= crej <= hi and cks > ALPHA_TEST) and (a.method == "pearson" or n >= 30):
                    refused[name] = True
    except Exception as e:
        print("gate_corr_calibration: harness failure: %s" % e)
        return 2
    for name, was in refused.items():
        required = a.method == "pearson" or name == "one-sided"
        lines.append("control %s: %s%s" % (name, "refused by the calibration test" if was else "NOT refused", "" if required else " (reported; at the gated sizes of Spearman only the one-sided tail is far enough from Student's t)"))
        if not was and required:
            bad.append("control %s was not refused: the gate cannot tell it from the program" % name)
    # C-F: power at rho = 0.5, n = 50 (Pearson only)
    if a.method != "pearson":
        print("gate_corr_calibration (%s, %d null replicates per n, seed %d):" % (a.method, a.replicates, a.seed))
        for ln in lines:
            print("  " + ln)
        print("gate_corr_calibration: %d problems" % len(bad))
        for b in bad:
            print("  " + b)
        return 1 if bad else 0
    n, rho = 50, 0.5
    tc = ss.t.ppf(0.975, n - 2)
    rc = tc / math.sqrt(tc * tc + n - 2)
    total = mpmath.quad(lambda x: fisher_density(x, rho, n), [-1, 0, rc, 1])
    power = mpmath.quad(lambda x: fisher_density(x, rho, n), [rc, 1]) + mpmath.quad(lambda x: fisher_density(x, rho, n), [-1, -rc])
    if abs(total - 1) > mpmath.mpf("1e-12"):
        bad.append("Fisher's density integrates to %s, not 1" % mpmath.nstr(total, 15))
    pr = SM64(a.seed * 7919)
    reps = a.power_replicates
    jobs = [sample(pr, n, rho) for _ in range(reps)]
    try:
        res = run_many(a.bin, jobs, a.workers)
    except Exception as e:
        print("gate_corr_calibration: harness failure: %s" % e)
        return 2
    rej = sum(1 for _, p in res if p < 0.05)
    plo = int(ss.binom.ppf(ALPHA_TEST / 2, reps, float(power)))
    phi = int(ss.binom.ppf(1 - ALPHA_TEST / 2, reps, float(power)))
    ok = plo <= rej <= phi
    lines.append("power rho=0.5 n=50: exact %.6f, rejected %d of %d (interval %d..%d)%s" % (float(power), rej, reps, plo, phi, "" if ok else "  <-- FAILS"))
    if not ok:
        bad.append("power")
    print("gate_corr_calibration (%d null replicates per n, seed %d):" % (a.replicates, a.seed))
    for ln in lines:
        print("  " + ln)
    print("gate_corr_calibration: %d problems" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
