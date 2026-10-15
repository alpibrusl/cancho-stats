#!/usr/bin/env python3
"""G5 / section 7.3 of docs/numerics.md: mutants of the special functions. A mutant is (name, file, old, new, ranges): `old` occurs exactly once in `file` and differs from `new`. Each is applied
to a COPY of lib/ and tools/ (with the installed stores), the driver is rebuilt, and the fast gates are rerun on the ranges the mutated function feeds: the accuracy sweeps (gate_special.py --fast
on the cached references), the domain and edge gate (gate_domain.py), and the thread check. Killed = a gate fails, the driver traps, or it hangs. A mutant that does not compile is INVALID, not killed.

    python3 scripts/mutants_special.py [--check] [--only SUBSTRING] [--jobs N]

`--check` only verifies that every pattern applies once (no compiler). The controls of gate G2 are the mutants whose name starts with `control-`: small, deliberate perturbations of one coefficient
or branch, which the sweeps must reject. A survivor is listed in EXPLAINED with the reason it cannot be killed (an equivalent mutant: it changes no result on any input); a survivor that is not
explained is a hole in the sweeps, and the sweep is extended, not the survivor explained away. Exit 1 if a mutant survives unexplained or is invalid.
"""
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
PY = sys.executable
REFS = str(ROOT / "build" / "refs-1-fast.json")

LG = "lib/special/lgam.cho"
IG = "lib/special/igam.cho"
IB = "lib/special/ibeta.cho"
DI = "lib/special/dist.cho"
XP = "lib/special/xprod.cho"
DD = "lib/special/dd.cho"
PR = "tools/special-check/parrun.cho"

# ranges (names or function names of gate_special.RANGES) each mutant is judged on
R_LG = "lgamma,lgamma1p,gamma_p,gamma_q,beta_lo,beta_up"
R_ERF = "erf,erfc,norm_cdf,norm_sf"
R_GAM = "gamma_p,gamma_q,chi2_cdf,chi2_sf"
R_BETA = "beta_lo,beta_up,t_cdf,t_sf,f_cdf,f_sf,binom_cdf,binom_sf"
R_Q = "qnorm"

MUTANTS = [
    # lgamma, delta, lgamma1p
    ("control-stirling-coefficient", LG, "let s = 1.0 / 12.0 - w", "let s = 1.0 / 12.001 - w", R_LG),
    ("lgamma-sign-of-x", LG, "return (x - 0.5) * math.log(x) - x + 0.9189385332046727 + delta(x);", "return (x - 0.5) * math.log(x) + x + 0.9189385332046727 + delta(x);", "lgamma"),
    ("lgamma-stirling-switch-too-low", LG, "    if x >= 10.0 {\n        return (x - 0.5)", "    if x >= 3.0 {\n        return (x - 0.5)", "lgamma"),
    ("lgamma-recurrence-product-sum", LG, "p = p * z;", "p = p + z;", "lgamma"),
    ("lgamma-recurrence-stops-early", LG, "while z < 10.0 {", "while z < 4.0 {", "lgamma"),
    ("lgamma-log-sqrt-2pi-truncated", LG, "return (z - 0.5) * math.log(z) - z + 0.9189385332046727 + delta(z) - math.log(p);", "return (z - 0.5) * math.log(z) - z + 0.918938533204 + delta(z) - math.log(p);", "lgamma"),
    ("lgamma-delta-drops-a-term", LG, "w * (691.0 / 360360.0 - w * (1.0 / 156.0))", "w * (691.0 / 360360.0)", "lgamma"),
    ("lgamma1p-euler-constant-truncated", LG, "0.0 - 0.5772156649015329 * a", "0.0 - 0.577215664901 * a", "lgamma1p"),
    ("lgamma1p-log1p-term-sign", LG, "(a - math.log1p(a)) + s", "(math.log1p(a) - a) + s", "lgamma1p"),
    ("lgamma1p-alternating-sign-dropped", LG, "sign = 0.0 - sign;", "sign = sign;", "lgamma1p"),
    ("lgamma1p-series-cut-short", LG, "    while k < 18 {\n        s = s + sign", "    while k < 5 {\n        s = s + sign", "lgamma1p"),
    # stirlerr, bd0, log_ratio
    ("stirlerr-series-coefficient", LG, "(1.0 / 12.0 - (1.0 / 360.0 -", "(1.0 / 12.0 - (1.0 / 361.0 -", "gamma_p,gamma_q,beta_lo,beta_up"),
    ("stirlerr-switch-too-low", LG, "if n <= 15.0 {", "if n <= 5.0 {", "gamma_p,gamma_q,chi2_cdf"),
    ("log-ratio-correction-dropped", LG, "return math.log(r) + e / x;", "return math.log(r);", "gamma_p,gamma_q,beta_lo,beta_up"),
    ("bd0-series-switch-window", LG, "if math.fabs(d) < 0.1 * (x + nph) {", "if math.fabs(d) < 0.9 * (x + nph) {", "gamma_p,gamma_q"),
    ("bd0-exact-np-low-part-dropped", LG, "return x * (log_ratio(x, nph) - npl / nph) + (nph - x) + npl;", "return x * log_ratio(x, nph) + (nph - x);", "beta_lo,beta_up"),
    # erf and normal
    ("exp-neg-sq-split-dropped", DI, "    let xl = x - xh;\n", "    let xl = 0.0;\n", R_ERF),
    ("exp-neg-sq-cross-term", DI, "math.exp(0.0 - k * xl * (x + xh))", "math.exp(0.0 - k * xl * x)", R_ERF),
    ("erf-series-switch-too-high", DI, "    let pr = exp_neg_sq(w, 1.0) * 2.0 * w / 1.7724538509055159;\n    if s < 1.5 {", "    let pr = exp_neg_sq(w, 1.0) * 2.0 * w / 1.7724538509055159;\n    if s < 30.0 {", "erfc"),
    ("erf-sqrt-pi-truncated", DI, "2.0 * w / 1.7724538509055159", "2.0 * w / 1.77245385", "erf,erfc"),
    ("erf-underflow-switch-too-low", DI, "    if w > 27.5 {", "    if w > 20.0 {", "erfc"),
    ("erfc-negative-branch-sign", DI, "        return (1.0 + p, 0);", "        return (1.0 - p, 0);", "erfc"),
    ("erf-odd-symmetry-dropped", DI, "    if x < 0.0 {\n        return (0.0 - p, 0);\n    }\n    return (p, 0);", "    return (p, 0);", "erf"),
    ("norm-tail-factor", DI, "let tail = 0.25 * pr * h;", "let tail = 0.5 * pr * h;", "norm_cdf,norm_sf"),
    ("norm-sqrt-two-over-pi-truncated", DI, "* w * 0.7978845608028654;", "* w * 0.79788456;", "norm_cdf,norm_sf"),
    ("norm-underflow-switch-too-low", DI, "    if w > 40.0 {", "    if w > 30.0 {", "norm_sf,norm_cdf"),
    ("norm-tails-swapped", DI, "    if z < 0.0 {\n        return (tail, big, 0);\n    }\n    return (big, tail, 0);", "    if z < 0.0 {\n        return (big, tail, 0);\n    }\n    return (tail, big, 0);", "norm_cdf,norm_sf"),
    ("norm-series-argument", DI, "    let s = w * w / 2.0;\n    let pr = exp_neg_sq(w, 0.5)", "    let s = w * w / 2.1;\n    let pr = exp_neg_sq(w, 0.5)", "norm_cdf,norm_sf"),
    # quantile
    ("control-no-newton-steps", DI, "    while i < 3 {\n        let lphi", "    while i < 0 {\n        let lphi", R_Q),
    ("qnorm-acklam-constant-and-no-steps", DI, "-3.969683028665376e+01 * r", "-3.96e+01 * r", R_Q),
    ("qnorm-tail-switch", DI, "    if p < 0.02425 {", "    if p < 0.0 {", R_Q),
    ("qnorm-upper-half-sign", DI, "        return (0.0 - qnorm_lower(1.0 - p), 0);", "        return (qnorm_lower(1.0 - p), 0);", R_Q),
    ("qnorm-newton-sign", DI, "        x = x - f / fp;", "        x = x + f / fp;", R_Q),
    ("qnorm-log-constant", DI, "- s + math.log(h) - 1.612085713764618;", "- s + math.log(h) - 1.6120857;", R_Q),
    ("qnorm-log-domain-switch", DI, "    if w > 5.0 {\n        let s = w * w / 2.0;", "    if w > 50.0 {\n        let s = w * w / 2.0;", R_Q),
    # incomplete gamma
    ("gamma-loader-switch-too-high", IG, "    if a >= 15.0 {\n        return math.exp(0.0 - lgam.stirlerr(a)", "    if a >= 1.0e6 {\n        return math.exp(0.0 - lgam.stirlerr(a)", "gamma_p,gamma_q"),
    ("gamma-series-tolerance", IG, "if term < sum * 1.0e-17 {", "if term < sum * 1.0e-5 {", R_GAM),
    ("gamma-cf-tolerance", IG, "if math.fabs(del - 1.0) < 2.5e-16 {", "if math.fabs(del - 1.0) < 2.5e-6 {", R_GAM),
    ("gamma-cf-term-sign", IG, "let an = 0.0 - fi * (fi - a);", "let an = fi * (fi - a);", R_GAM),
    ("gamma-cf-start-form", IG, "var b = (x - a) + 1.0;", "var b = x + 1.0 - a + 0.0 * x;", "gamma_p,gamma_q"),
    ("gamma-switch-point", IG, "    if a >= 0.1 && x < a + 1.0 {", "    if a >= 0.1 && x < 2.0 * a + 1.0 {", "gamma_p,gamma_q"),
    ("gamma-small-a-path-disabled", IG, "    if a < 0.1 && x < 1.0 {", "    if a < 0.0 && x < 1.0 {", "gamma_p,gamma_q"),
    ("gamma-small-a-expm1", IG, "clamp01(0.0 - math.expm1(l) - el * a * sum)", "clamp01(1.0 - el - el * a * sum)", "gamma_q"),
    ("gamma-upper-from-prefactor-missing-a", IG, "let q = clamp01(pr * a * h);", "let q = clamp01(pr * h);", R_GAM),
    ("gamma-p-from-q", IG, "return (clamp01(1.0 - q), q, 0);", "return (q, q, 0);", "gamma_p"),
    ("gamma-prefactor-underflow-branch", IG, "        return (1.0, 0.0, 0);\n    }\n    if a >= 0.1 && x < a + 1.0 {", "        return (0.0, 1.0, 0);\n    }\n    if a >= 0.1 && x < a + 1.0 {", "gamma_p,gamma_q"),
    ("gamma-cap-one", IG, "    return 2000 + truncate(100.0 * sqrt(a));", "    return 1;", "gamma_p,gamma_q"),
    ("gamma-prefactor-direct-form", IG, "return math.exp(a * math.log(x) - x - lgam.lgamma1p(a));", "return math.exp(a * math.log(x) - x);", "gamma_p,gamma_q"),
    # incomplete beta
    ("beta-betaln-switch", IB, "    if q < 10.0 {\n        return lgam.lgamma_pos(p)", "    if q < 1.0 {\n        return lgam.lgamma_pos(p)", "beta_lo,beta_up,t_cdf"),
    ("beta-betaln-bracket-sign", IB, "- p * math.log(s) + p +", "+ p * math.log(s) + p +", "beta_lo,beta_up,t_cdf"),
    ("beta-front-log-branches-swapped", IB, "        lx = math.log(xh);\n        ly = math.log1p(0.0 - xh);", "        lx = math.log1p(0.0 - xh);\n        ly = math.log(xh);", "beta_lo,beta_up"),
    ("beta-cf-tolerance", IB, "if math.fabs((delh - 1.0) + dell) < 2.5e-16 {", "if math.fabs((delh - 1.0) + dell) < 2.5e-6 {", R_BETA),
    ("beta-cf-first-term-uses-x", IB, "let (ayh, ayl) = dd.mul(a, 0.0, yh, yl);", "let (ayh, ayl) = dd.mul(a, 0.0, xh, xl);", R_BETA),
    ("beta-cf-even-term", IB, "let (e1h, e1l) = dd.two_sum(a, m2 - 1.0);", "let (e1h, e1l) = dd.two_sum(a, m2);", R_BETA),
    ("beta-cf-odd-term-sign", IB, "lentz_step(0.0 - bbh, 0.0 - bbl, sd1h, sd1l, sc1h, sc1l);", "lentz_step(bbh, bbl, sd1h, sd1l, sc1h, sc1l);", R_BETA),
    ("beta-lentz-one-dropped", IB, "let (d1h, d1l) = dd.add(1.0, 0.0, t1h, t1l);", "let (d1h, d1l) = dd.add(0.0, 0.0, t1h, t1l);", R_BETA),
    ("beta-switch-point", IB, "    if xh < (a + 1.0) / (a + b + 2.0) {", "    if xh < 0.5 {", R_BETA),
    ("beta-complement-recompute-off", IB, "        if up < 0.3 && up * sqrt(maxab(a, b)) < 1.0 {", "        if up < 0.0 {", "beta_lo,beta_up,binom_cdf,binom_sf"),
    ("beta-lower-from-upper-dropped", IB, "        lo = 1.0 - up;\n        if lo < 0.3", "        lo = up;\n        if lo < 0.3", "beta_lo,beta_up"),
    ("beta-cap-one", IB, "    return 2000 + truncate(100.0 * sqrt(m));", "    return 1;", "beta_lo,beta_up"),
    ("beta-loader-switch", IB, "    if a >= 10.0 && b >= 10.0 {\n        let n = a + b;", "    if a >= 10.0 || b >= 10.0 {\n        let n = a + b;", "beta_lo,beta_up"),
    ("beta-exact-complement-dropped", IB, "        let (h1, l1) = dd.two_sum(1.0, 0.0 - x);\n        yh = h1;\n        yl = l1;", "        yh = 1.0 - x;", "t_cdf,t_sf,beta_lo,beta_up"),
    # double-double and the two-product
    ("dd-add-low-parts-dropped", DD, "    let (t1, t2) = two_sum(al, bl);\n    let (u1, u2) = quick_two_sum(s1, s2 + t1);\n    return quick_two_sum(u1, u2 + t2);", "    let (u1, u2) = quick_two_sum(s1, s2);\n    return (u1, u2);", "beta_lo,beta_up,t_cdf"),
    ("dd-mul-cross-terms-dropped", DD, "    return quick_two_sum(p, e + (ah * bl + al * bh));", "    return quick_two_sum(p, e);", "beta_lo,beta_up,t_cdf"),
    ("dd-div-correction-dropped", DD, "    let q3 = r2h / bh;", "    let q3 = 0.0;", "beta_lo,beta_up,t_cdf"),
    ("xprod-split-constant", XP, "let c = 134217729.0 * a;", "let c = 134217730.0 * a;", "gamma_p,gamma_q,beta_lo,beta_up,t_cdf"),
    ("xprod-error-term", XP, "return (p, ah * bh - p + ah * bl + al * bh + al * bl);", "return (p, ah * bh - p + ah * bl + al * bh);", "gamma_p,gamma_q,beta_lo,beta_up,t_cdf"),
    # the distributions
    ("t-tail-halving-dropped", DI, "    let tail = lo / 2.0;", "    let tail = lo;", "t_cdf,t_sf"),
    ("t-other-tail", DI, "    let other = (1.0 + up) / 2.0;", "    let other = up / 2.0;", "t_cdf,t_sf"),
    ("t-sign-branch", DI, "    if t > 0.0 {\n        return (other, tail, 0);\n    }\n    return (tail, other, 0);", "    if t > 0.0 {\n        return (tail, other, 0);\n    }\n    return (other, tail, 0);", "t_cdf,t_sf"),
    ("t-argument", DI, "        x = nu / (nu + at * at);", "        x = nu / (nu + at);", "t_cdf,t_sf"),
    ("t-large-argument-branch", DI, "        let r = nu / at;\n        x = r / (r + at);\n        y = at / (r + at);", "        let r = nu / at;\n        x = r / (r + at);\n        y = 1.0 - x;", "t_cdf,t_sf"),
    ("chi2-half-argument", DI, "igam.gamma_pq(nu / 2.0, x / 2.0)", "igam.gamma_pq(nu / 2.0, x)", "chi2_cdf,chi2_sf"),
    ("f-argument", DI, "        x = r / (r + d2);\n        y = d2 / (r + d2);", "        x = r / (r + d1);\n        y = d2 / (r + d2);", "f_cdf,f_sf"),
    ("f-parameters-swapped", DI, "ibeta.beta_xy(d1 / 2.0, d2 / 2.0, x, y)", "ibeta.beta_xy(d2 / 2.0, d1 / 2.0, x, y)", "f_cdf,f_sf"),
    ("binom-parameters-swapped", DI, "ibeta.beta_xy(n - k, k + 1.0, q, p)", "ibeta.beta_xy(k + 1.0, n - k, q, p)", "binom_cdf,binom_sf"),
    ("binom-q", DI, "    var q = 1.0 - p;", "    var q = p;", "binom_cdf,binom_sf"),
    ("binom-k-offset", DI, "ibeta.beta_xy(n - k, k + 1.0, q, p)", "ibeta.beta_xy(n - k, k, q, p)", "binom_cdf,binom_sf"),
    # domains: the refusals (the domain gate)
    ("lgamma-domain-zero", DI, "    if x <= 0.0 {\n        return (nan(), rules.domain_param());\n    }\n    if x > 1.0e300", "    if x < 0.0 {\n        return (nan(), rules.domain_param());\n    }\n    if x > 1.0e300", "lgamma"),
    ("qnorm-domain-edges", DI, "    if p <= 0.0 || p >= 1.0 {", "    if p < 0.0 || p > 1.0 {", "qnorm"),
    ("refusal-tag-of-param", DI, "    if nu <= 0.0 || is_inf(nu) {\n        return (nan(), nan(), rules.domain_param());", "    if nu <= 0.0 || is_inf(nu) {\n        return (nan(), nan(), rules.domain_prob());", "t_cdf"),
    ("gamma-range-limit", IG, "    if a > 1.0e8 {", "    if a > 1.0e12 {", "gamma_p"),
    # threads
    ("threads-block-slot", PR, "copy_ints(ov[lo..hi], contents(ow), hi - lo);", "copy_ints(ov[0..hi - lo], contents(ow), hi - lo);", "erf"),
    ("threads-partition", PR, "    let lo = i * n / t;", "    let lo = i * n / (t + 1);", "erf"),
]

# a survivor that no input can kill: the mutation changes no result. Name -> reason.
EXPLAINED = {}


def problems():
    out = []
    names = set()
    cache = {}
    for name, file, old, new, ranges in MUTANTS:
        if name in names:
            out.append("mutant %s: named twice" % name)
        names.add(name)
        p = ROOT / file
        if not p.exists():
            out.append("mutant %s: no such file %s" % (name, file))
            continue
        t = cache.setdefault(file, p.read_text())
        n = t.count(old)
        if n != 1:
            out.append("mutant %s: pattern occurs %d times in %s (must be exactly 1)" % (name, n, file))
        if old == new:
            out.append("mutant %s: replacement equals pattern" % name)
    return out


def checks(driver, tmp, ranges):
    """-> (killed, reason)"""
    env = dict(os.environ)
    gs = subprocess.run([PY, str(ROOT / "scripts" / "gate_special.py"), "--driver", driver, "--fast", "--ranges", ranges, "--refs", REFS, "--workers", "1"], capture_output=True, text=True, timeout=1800, env=env)
    if gs.returncode != 0:
        return True, "accuracy sweep (exit %d)" % gs.returncode
    gd = subprocess.run([PY, str(ROOT / "scripts" / "gate_domain.py"), "--driver", driver, "--fuzz", "1500", "--threads", "4"], capture_output=True, text=True, timeout=1800, env=env)
    if gd.returncode != 0:
        return True, "domain and fuzz gate (exit %d)" % gd.returncode
    # thread check: the bytes of --threads T equal the sequential bytes
    sys.path.insert(0, str(ROOT / "scripts"))
    import gate_special as g

    cs = [(g.SM64(5).uni(-3, 3) + i * 0.01,) for i in range(300)]
    base = g.run_driver(driver, "erf", cs)
    for t in (2, 3, 7):
        if g.run_driver(driver, "erf", cs, threads=t) != base:
            return True, "thread check (--threads %d differs)" % t
    return False, ""


def main():
    bad = problems()
    if bad:
        print("\n".join(bad))
        return 1
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else ""
    if "--check" in sys.argv:
        print("mutants_special: %d mutants, every pattern applies once" % len(MUTANTS))
        return 0
    if not pathlib.Path(REFS).exists():
        print("mutants_special: no cached references %s (run scripts/gate_special.py --fast first)" % REFS)
        return 2
    survivors, invalid, killed, explained = [], [], 0, 0
    for name, file, old, new, ranges in MUTANTS:
        if only and only not in name:
            continue
        with tempfile.TemporaryDirectory() as t:
            tmp = pathlib.Path(t)
            for d in ("lib", "tools"):
                shutil.copytree(ROOT / d, tmp / d)
            shutil.copy(ROOT / "cancho.toml", tmp / "cancho.toml")
            (tmp / "build").mkdir()
            shutil.copytree(ROOT / "build" / "deps", tmp / "build" / "deps")
            f = tmp / file
            f.write_text(f.read_text().replace(old, new, 1))
            r = subprocess.run([CANCHO, "build", "--bin", "special-check"], cwd=tmp, capture_output=True, text=True, timeout=600)
            if r.returncode != 0:
                invalid.append(name)
                print("%-44s INVALID (does not compile): %s" % (name, (r.stderr or r.stdout)[-200:].replace("\n", " ")), flush=True)
                continue
            try:
                dead, why = checks(str(tmp / "build" / "special-check"), tmp, ranges)
            except subprocess.TimeoutExpired:
                dead, why = True, "hang (timeout)"
            except Exception as e:  # a driver that traps makes a gate raise: that is a kill
                dead, why = True, "trap or crash (%s)" % str(e)[:60]
        if dead:
            killed += 1
            print("%-44s killed by %s" % (name, why), flush=True)
        elif name in EXPLAINED:
            explained += 1
            print("%-44s SURVIVED (explained: %s)" % (name, EXPLAINED[name]), flush=True)
        else:
            survivors.append(name)
            print("%-44s SURVIVED" % name, flush=True)
    print("mutants_special: %d killed, %d survived and explained, %d unexplained, %d invalid" % (killed, explained, len(survivors), len(invalid)))
    if survivors:
        print("SURVIVORS: " + ", ".join(survivors))
    return 1 if (survivors or invalid) else 0


if __name__ == "__main__":
    sys.exit(main())
