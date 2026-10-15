#!/usr/bin/env python3
"""G1 (docs/numerics.md section 7.1): the special functions against `mpmath` at high precision, over the ranges of section 3.2, with the bounds of that table.

    python3 scripts/gate_special.py [--driver build/special-check] [--n-scale 1.0] [--fast] [--ranges NAME,NAME] [--workers 8] [--seed 1] [--refs FILE]

For every range: draw the cases from a seeded SplitMix64 stream (the inputs are exact doubles), compute the reference at high precision (scripts/special_ref.py), run the driver on the
cases (stdin: `id hex hex ...`, stdout: `id hex code`), and check each result against its bound B(f) = B0 + 4e-16 |ln f| (relative), B0 max(1, |f|) (lgamma, absolute) or 1e-15 + 1e-13 |x| (qnorm).
A true value below the smallest normal double is judged by an absolute bound of 2e-323; a result that is NaN, infinite, or accompanied by a refusal code is a failure inside a range.

The references are cached (`--refs`, default build/refs-<seed>-<scale>.json) so that the mutants, which rerun the driver only, do not recompute them.
Exit 0: every case of every range is inside its bound. Exit 1: a case is outside (the report names the worst case per range). Exit 2: the harness failed.
"""
import argparse
import hashlib
import json
import math
import multiprocessing
import os
import pathlib
import struct
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

M64 = (1 << 64) - 1


class SM64:
    """SplitMix64, the generator of cancho-ml's docs/rng.md, so a case list is a function of the seed alone."""

    def __init__(self, seed):
        self.s = seed & M64

    def next(self):
        self.s = (self.s + 0x9E3779B97F4A7C15) & M64
        z = self.s
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
        return z ^ (z >> 31)

    def u(self):
        return (self.next() >> 11) / float(1 << 53)

    def uni(self, lo, hi):
        return lo + (hi - lo) * self.u()

    def logu(self, lo, hi):
        return math.exp(math.log(lo) + self.u() * (math.log(hi) - math.log(lo)))

    def pick(self, seq):
        return seq[self.next() % len(seq)]


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def hx(x):
    return "%016x" % bits(x)


def unhx(h):
    return struct.unpack("<d", struct.pack("<Q", int(h, 16)))[0]


def nb(x, k):
    """k steps of nextafter from x (k may be negative)."""
    for _ in range(abs(k)):
        x = math.nextafter(x, math.inf if k > 0 else -math.inf)
    return x


def around(x):
    return [nb(x, -1), x, nb(x, 1)]


# ------------------------------------------------------------------------------------------------------------- the case generators
# gen(rng, i) -> tuple of doubles; i is the index of the case, its low bits choose the mode so every mode is represented in any prefix.


def g_lgamma(r, i):
    m = i % 5
    if m == 0:
        return (r.logu(1e-300, 1e15),)
    if m == 1:
        return (r.uni(0.5, 12.0),)
    if m == 2:
        return (r.pick([1.0, 2.0]) + (r.u() - 0.5) * 1e-3,)
    if m == 3:
        return (r.pick(around(10.0)),)
    return (r.pick([nb(10.0, -3), nb(10.0, 3), 9.99999, 10.00001, 0.5, 1.5, 3.0, 20.0]) * (1 + 0 * r.u()),)


def g_erf(r, i):
    m = i % 3
    s = 1.0 if r.u() < 0.5 else -1.0
    if m == 0:
        return (s * r.logu(1e-300, 6.0),)
    if m == 1:
        return (r.uni(-6.0, 6.0),)
    return (s * r.pick(around(1e-8) + around(1.0) + around(0.5)),)


def g_erfc(r, i):
    m = i % 3
    if m == 0:
        return (r.uni(-6.0, 6.0),)
    if m == 1:
        return (r.uni(6.0, 26.5),)
    return (r.pick(around(0.5) + around(1.0) + around(2.0) + around(26.5) + [26.0, 10.0, -6.0]),)


def g_gamma_mid(r, i):
    a = r.logu(1e-2, 1e4)
    m = i % 5
    if m == 0:
        x = r.logu(1e-300, 1e5)
    elif m == 1 or m == 2:
        x = a + r.uni(-12.0, 12.0) * math.sqrt(a)
        if x <= 0:
            x = a * r.u()
    elif m == 3:
        a = r.pick([15.0, nb(15.0, -1), nb(15.0, 1), 14.999, 15.001])
        x = a * r.uni(0.5, 1.5)
    else:
        x = r.pick(around(a + 1))
    return (a, max(x, 1e-300))


def g_gamma_big(r, i):
    a = r.logu(1e4, 1e7)
    m = i % 4
    if m == 3:
        x = r.pick(around(a + 1))
    else:
        x = a + r.uni(-12.0, 12.0) * math.sqrt(a)
    return (a, x)


def g_gamma_tiny(r, i):
    return (r.logu(1e-300, 1e-2), r.logu(1e-300, 50.0))


def _beta_mean_sd(a, b):
    s = a + b
    return a / s, math.sqrt(a * b / (s * s * (s + 1)))


def g_beta_small(r, i):
    a = r.logu(1e-2, 1e3)
    b = r.logu(1e-2, 1e3)
    m = i % 7
    if m == 0:
        x = r.u()
    elif m == 1 or m == 2:
        mu, sd = _beta_mean_sd(a, b)
        x = mu + r.uni(-12.0, 12.0) * sd
    elif m == 3:
        x = r.logu(1e-300, 1e-3)
    elif m == 4:
        x = 1.0 - r.logu(1e-16, 1e-3)
    elif m == 5:
        x = r.pick(around((a + 1) / (a + b + 2)))
    else:
        a = r.pick(around(10.0))
        b = r.pick(around(10.0) + [20.0])
        mu, sd = _beta_mean_sd(a, b)
        x = mu + r.uni(-8.0, 8.0) * sd
    if not (0.0 < x < 1.0):
        x = r.uni(0.01, 0.99)
    return (a, b, x)


def g_beta_large(r, i):
    big = r.logu(1e3, 1e6)
    other = r.logu(1e-2, 1e6)
    a, b = (big, other) if (i & 1) == 0 else (other, big)
    mu, sd = _beta_mean_sd(a, b)
    if i % 5 == 4:
        x = r.pick(around((a + 1) / (a + b + 2)))
    else:
        x = mu + r.uni(-12.0, 12.0) * sd
    if not (0.0 < x < 1.0):
        x = min(max(mu, 1e-300), 1.0 - 2 ** -53)
    return (a, b, x)


def g_norm(r, i):
    m = i % 3
    if m == 0:
        return (r.uni(-38.0, 38.0),)
    if m == 1:
        return (r.uni(-1.0, 1.0),)
    return (r.uni(-38.0, -5.0) * (1.0 if r.u() < 0.5 else -1.0),)


def g_qnorm(r, i):
    m = i % 5
    if m == 0:
        return (r.logu(5e-324, 0.5),)
    if m == 1:
        return (r.u() or 0.25,)
    if m == 2:
        return (nb(0.5, r.next() % 41 - 20),)
    if m == 3:
        return (r.pick(around(0.02425) + around(0.97575)),)
    return (r.logu(1e-300, 1e-3) if r.u() < 0.5 else 1.0 - r.logu(1e-16, 1e-2),)


def g_t(r, i):
    nu = r.logu(0.1, 1e6)
    if i % 10 == 9:
        nu = 1e8
    tmax = 60.0 if nu > 1e3 else 1e4
    t = r.logu(1e-6, tmax) * (1.0 if r.u() < 0.5 else -1.0)
    if i % 11 == 10:
        t = 0.0
    return (t, nu)


def g_chi2_mid(r, i):
    nu = r.logu(1e-2, 1e4)
    if i % 2 == 0:
        x = max(nu + r.uni(-12.0, 12.0) * math.sqrt(2 * nu), 1e-300)
    else:
        x = r.logu(1e-300, 5 * nu + 50)
    return (x, nu)


def g_chi2_big(r, i):
    nu = r.logu(1e4, 1e6)
    return (nu + r.uni(-12.0, 12.0) * math.sqrt(2 * nu), nu)


def g_f(r, i):
    d1 = r.logu(0.1, 1e5)
    d2 = r.logu(0.1, 1e5)
    if i % 3 == 0:
        f = r.logu(1e-3, 1e3)
    else:
        # near the centre of the distribution: x = d1 f / (d1 f + d2) near its mean
        mu, sd = _beta_mean_sd(d1 / 2, d2 / 2)
        x = mu + r.uni(-12.0, 12.0) * sd
        if not (0.0 < x < 1.0):
            x = mu
        f = d2 * x / (d1 * (1 - x))
    return (f, d1, d2)


def g_binom(r, i):
    n = float(max(1, int(round(r.logu(1.0, 1e6)))))
    p = r.logu(1e-6, 0.5)
    if i & 1:
        p = 1.0 - p
    mu = n * p
    sd = math.sqrt(max(n * p * (1 - p), 1e-12))
    k = round(mu + r.uni(-8.0, 8.0) * sd) if i % 3 else int(r.u() * n)
    k = float(min(max(k, 0), n - 1))
    return (k, n, p)


# (name, fn, n cases, generator, kind, B0, slow)
RANGES = [
    ("lgamma", "lgamma", 20000, g_lgamma, "lgamma", 1e-14, False),
    ("erf", "erf", 20000, g_erf, "rel", 2e-13, False),
    ("erfc", "erfc", 20000, g_erfc, "rel", 2e-13, False),
    ("gamma_p mid", "gamma_p", 20000, g_gamma_mid, "rel", 2e-13, False),
    ("gamma_q mid", "gamma_q", 20000, g_gamma_mid, "rel", 2e-13, False),
    ("gamma_p big", "gamma_p", 2000, g_gamma_big, "rel", 5e-12, True),
    ("gamma_q big", "gamma_q", 2000, g_gamma_big, "rel", 5e-12, True),
    ("gamma_p tiny", "gamma_p", 20000, g_gamma_tiny, "rel", 2e-12, False),
    ("gamma_q tiny", "gamma_q", 20000, g_gamma_tiny, "rel", 2e-12, False),
    ("beta_lo small", "beta_lo", 20000, g_beta_small, "rel", 5e-13, False),
    ("beta_up small", "beta_up", 20000, g_beta_small, "rel", 5e-13, False),
    ("beta_lo large", "beta_lo", 2000, g_beta_large, "rel", 1e-11, True),
    ("beta_up large", "beta_up", 2000, g_beta_large, "rel", 1e-11, True),
    ("norm_cdf", "norm_cdf", 20000, g_norm, "rel", 2e-13, False),
    ("norm_sf", "norm_sf", 20000, g_norm, "rel", 2e-13, False),
    ("qnorm", "qnorm", 20000, g_qnorm, "qnorm", 0.0, False),
    ("t_cdf", "t_cdf", 20000, g_t, "rel", 1e-12, False),
    ("t_sf", "t_sf", 20000, g_t, "rel", 1e-12, False),
    ("chi2_cdf mid", "chi2_cdf", 20000, g_chi2_mid, "rel", 2e-13, False),
    ("chi2_sf mid", "chi2_sf", 20000, g_chi2_mid, "rel", 2e-13, False),
    ("chi2_cdf big", "chi2_cdf", 2000, g_chi2_big, "rel", 2e-12, True),
    ("chi2_sf big", "chi2_sf", 2000, g_chi2_big, "rel", 2e-12, True),
    ("f_cdf", "f_cdf", 20000, g_f, "rel", 1e-12, False),
    ("f_sf", "f_sf", 20000, g_f, "rel", 1e-12, False),
    ("binom_cdf", "binom_cdf", 20000, g_binom, "rel", 5e-12, False),
    ("binom_sf", "binom_sf", 20000, g_binom, "rel", 5e-12, False),
]

B1 = 4e-16
SMALLEST_NORMAL = 2.2250738585072014e-308


def cases_of(rng_spec, seed, scale, fast):
    name, fn, n, gen, kind, b0, slow = rng_spec
    n = max(1, int(n * scale))
    if fast:
        n = 600 if not slow else 200
    h = int.from_bytes(hashlib.sha256(("%s/%d" % (name, seed)).encode()).digest()[:8], "little")
    r = SM64(h)
    return [gen(r, i) for i in range(n)]


def _ref(job):
    fn, args = job
    import special_ref
    from mpmath import mp

    try:
        v = special_ref.ref(fn, args)
        return mp.nstr(v, 45)
    except Exception as e:  # a reference that fails is reported, never guessed
        return "ERR:" + type(e).__name__ + ":" + str(e)[:60]


def references(specs, seed, scale, fast, workers, cache):
    jobs = []
    index = []
    for spec in specs:
        cs = cases_of(spec, seed, scale, fast)
        for k, c in enumerate(cs):
            jobs.append((spec[1], c))
            index.append((spec[0], k))
    store = {}
    if cache and pathlib.Path(cache).exists():
        store = json.load(open(cache))
    key = lambda j: j[0] + "|" + ",".join(hx(a) for a in j[1])
    todo = [j for j in jobs if key(j) not in store]
    if todo:
        t0 = time.time()
        with multiprocessing.Pool(workers) as pool:
            for j, v in zip(todo, pool.imap(_ref, todo, chunksize=16)):
                store[key(j)] = v
        print("references: %d computed in %.0f s (%d workers)" % (len(todo), time.time() - t0, workers), flush=True)
        if cache:
            pathlib.Path(cache).parent.mkdir(parents=True, exist_ok=True)
            json.dump(store, open(cache, "w"))
    return {i: store[key(j)] for i, j in zip(index, jobs)}, key


def run_driver(driver, fn, cases, timeout=900, threads=None):
    lines = []
    for k, c in enumerate(cases):
        lines.append("%d %s" % (k, " ".join(hx(a) for a in c)))
    argv = [driver, fn] + (["--threads", str(threads)] if threads else [])
    r = subprocess.run(argv, input="\n".join(lines) + "\n", capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError("driver %s exit %d: %s" % (fn, r.returncode, (r.stderr or r.stdout)[-300:]))
    out = [ln.split() for ln in r.stdout.strip().split("\n") if ln.strip()]
    if len(out) != len(cases):
        raise RuntimeError("driver %s answered %d lines for %d cases" % (fn, len(out), len(cases)))
    return out


def judge(kind, b0, y, code, tag, f):
    """(ok, ratio of error to bound, absolute error) for a double y against the reference f (mpf)."""
    from mpmath import mp, mpf

    if code != 0 or y != y or abs(y) == math.inf:
        return False, float("inf"), float("inf")
    mp.dps = 60
    err = abs(mpf(y) - f)
    af = abs(f)
    if kind == "lgamma":
        bound = mpf(b0) * max(mpf(1), af)
    elif kind == "qnorm":
        bound = mpf("1e-15") + mpf("1e-13") * af
    elif af < SMALLEST_NORMAL:
        bound = mpf("2e-323")
        if af == 0:
            return (y == 0.0), (0.0 if y == 0.0 else float("inf")), float(err)
    else:
        bound = (mpf(b0) + mpf(B1) * abs(mp.log(af))) * af
    ratio = float(err / bound)
    return err <= bound, ratio, float(err / af) if af != 0 else float(err)


def evaluate(driver, specs, seed, scale, fast, workers, cache, quiet=False):
    refs, key = references(specs, seed, scale, fast, workers, cache)
    from mpmath import mpf

    results = []
    for spec in specs:
        name, fn, n, gen, kind, b0, slow = spec
        cs = cases_of(spec, seed, scale, fast)
        out = run_driver(driver, fn, cs)
        worst = (-1.0, None)
        worst_rel = 0.0
        bad = 0
        first_bad = None
        errs = 0
        for k, (c, o) in enumerate(zip(cs, out)):
            rv = refs[(name, k)]
            if rv.startswith("ERR:"):
                errs += 1
                continue
            f = mpf(rv)
            y = unhx(o[1])
            ok, ratio, rel = judge(kind, b0, y, int(o[2]), None, f)
            if not ok:
                bad += 1
                if first_bad is None or ratio > first_bad[0]:
                    first_bad = (ratio, c, y, rv)
            if ratio > worst[0]:
                worst = (ratio, (c, y, rv))
            if rel != float("inf"):
                worst_rel = max(worst_rel, rel)
        results.append((name, len(cs), bad, worst, worst_rel, errs, first_bad))
        if not quiet:
            print("%-16s n=%-6d outside=%-5d max error/bound=%-9.3g max relative=%-9.3g reference-errors=%d%s" % (
                name, len(cs), bad, max(worst[0], 0.0), worst_rel, errs,
                "" if not bad else "  FIRST WORST: args=%s got=%r ref=%s" % (",".join(hx(a) for a in first_bad[1]), first_bad[2], first_bad[3][:25])), flush=True)
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", default=str(ROOT / "build" / "special-check"))
    ap.add_argument("--n-scale", type=float, default=1.0)
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--ranges", default="")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--refs", default="")
    a = ap.parse_args()
    specs = RANGES
    if a.ranges:
        want = set(a.ranges.split(","))
        specs = [s for s in RANGES if s[0] in want or s[1] in want]
        if not specs:
            print("no such range")
            return 2
    cache = a.refs or str(ROOT / "build" / ("refs-%d-%s.json" % (a.seed, "fast" if a.fast else ("%g" % a.n_scale))))
    try:
        res = evaluate(a.driver, specs, a.seed, a.n_scale, a.fast, a.workers, cache)
    except Exception as e:
        print("gate_special: harness failure: %s" % e)
        return 2
    bad = sum(r[2] for r in res)
    referr = sum(r[5] for r in res)
    print("gate_special: %d ranges, %d cases outside the bound, %d reference failures" % (len(res), bad, referr))
    if referr:
        return 2
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
