#!/usr/bin/env python3
"""G3 (docs/numerics.md section 7.1): an independent second reference. The sweeps of gate_special.py are judged against `mpmath` at high precision; this script puts `scipy.special` (double
precision, written by other people in another algorithm) beside it on the same cases, so that an error SHARED by the driver and the mpmath-based reference cannot hide.

    python3 scripts/gate_independence.py [--driver build/special-check] [--fast | --n-scale S] [--refs FILE] [--seed 1]

Per range it reports the largest relative disagreement of scipy with the mpmath reference and of the driver with scipy, and counts the points where the two REFERENCES disagree by more than
1e-10 relative. At such a point there are three cases, decided by who is isolated: (a) scipy is the outlier (the driver agrees with mpmath: scipy's own error, counted and listed by range,
never hidden); (b) the mpmath reference is the outlier (the driver agrees with scipy, not with mpmath: the gate FAILS, the reference is suspect); (c) the driver is the outlier: also covered by gate_special.py.
No sentence of this repository says that a result "matches" scipy: the numbers printed here are the disagreements, with the cases they were measured on.
Exit 0: no point of case (b). Exit 1: a point of case (b). Exit 2: the harness failed.
"""
import argparse
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import gate_special as g  # noqa: E402


def scipy_functions():
    from scipy import special as sp

    return {
        "lgamma": lambda a: sp.gammaln(a[0]),
        "erf": lambda a: sp.erf(a[0]),
        "erfc": lambda a: sp.erfc(a[0]),
        "gamma_p": lambda a: sp.gammainc(a[0], a[1]),
        "gamma_q": lambda a: sp.gammaincc(a[0], a[1]),
        "beta_lo": lambda a: sp.betainc(a[0], a[1], a[2]),
        "beta_up": lambda a: sp.betaincc(a[0], a[1], a[2]),
        "norm_cdf": lambda a: sp.ndtr(a[0]),
        "norm_sf": lambda a: sp.ndtr(-a[0]),
        "qnorm": lambda a: sp.ndtri(a[0]),
        "t_cdf": lambda a: sp.stdtr(a[1], a[0]),
        "t_sf": lambda a: sp.stdtr(a[1], -a[0]),
        "chi2_cdf": lambda a: sp.chdtr(a[1], a[0]),
        "chi2_sf": lambda a: sp.chdtrc(a[1], a[0]),
        "f_cdf": lambda a: sp.fdtr(a[1], a[2], a[0]),
        "f_sf": lambda a: sp.fdtrc(a[1], a[2], a[0]),
        "binom_cdf": lambda a: sp.bdtr(a[0], a[1], a[2]),
        "binom_sf": lambda a: sp.bdtrc(a[0], a[1], a[2]),
    }


def rel(x, y):
    if x == y:
        return 0.0
    if y == 0.0 or y != y or x != x:
        return math.inf
    return abs(x - y) / abs(y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", default=str(ROOT / "build" / "special-check"))
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--n-scale", type=float, default=1.0)
    ap.add_argument("--refs", default="")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    from mpmath import mpf

    sp = scipy_functions()
    specs = [s for s in g.RANGES if s[1] in sp]
    cache = a.refs or str(ROOT / "build" / ("refs-%d-%s.json" % (a.seed, "fast" if a.fast else ("%g" % a.n_scale))))
    refs, key = g.references(specs, a.seed, a.n_scale, a.fast, a.workers, cache)
    bad_b = []
    print("%-16s %6s | %-26s | %-26s | %s" % ("range", "n", "scipy vs mpmath: max rel", "driver vs scipy: max rel", "points > 1e-10 between the references: scipy isolated / mpmath isolated"))
    for spec in specs:
        name, fn, n, gen, kind, b0, slow = spec
        cs = g.cases_of(spec, a.seed, a.n_scale, a.fast)
        out = g.run_driver(a.driver, fn, cs)
        worst_sm = 0.0
        worst_ds = 0.0
        iso_s = 0
        iso_m = 0
        used = 0
        for k, (c, o) in enumerate(zip(cs, out)):
            rv = refs[(name, k)]
            if rv.startswith("ERR:") or int(o[2]) != 0:
                continue
            f = float(mpf(rv))
            y = g.unhx(o[1])
            try:
                s = float(sp[fn](c))
            except Exception:
                continue
            if s != s or math.isinf(s) or abs(f) < 1e-290 or abs(s) < 1e-290:
                continue  # below the range where a double-precision library is a reference at all
            used += 1
            r_sm = rel(s, f) if kind != "lgamma" else abs(s - f) / max(1.0, abs(f))
            r_ds = rel(y, s) if kind != "lgamma" else abs(y - s) / max(1.0, abs(s))
            if kind == "qnorm":
                r_sm = abs(s - f) / max(1e-3, abs(f))
                r_ds = abs(y - s) / max(1e-3, abs(s))
            worst_sm = max(worst_sm, r_sm)
            worst_ds = max(worst_ds, r_ds)
            if r_sm > 1e-10:
                r_dm = rel(y, f) if kind != "lgamma" else abs(y - f) / max(1.0, abs(f))
                if kind == "qnorm":
                    r_dm = abs(y - f) / max(1e-3, abs(f))
                if r_dm <= 1e-10:
                    iso_s += 1
                else:
                    iso_m += 1
                    bad_b.append((name, c, y, s, rv[:25]))
        print("%-16s %6d | %-26.3g | %-26.3g | %d / %d" % (name, used, worst_sm, worst_ds, iso_s, iso_m))
    for b in bad_b[:5]:
        print("  mpmath isolated: %s %s driver=%r scipy=%r mpmath=%s" % (b[0], ",".join(g.hx(x) for x in b[1]), b[2], b[3], b[4]))
    print("gate_independence: %s" % ("FAILED: %d points where the mpmath reference is isolated" % len(bad_b) if bad_b else "ok (no point where the mpmath reference is isolated)"))
    return 1 if bad_b else 0


if __name__ == "__main__":
    sys.exit(main())
