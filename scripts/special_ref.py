"""The reference values of the special-function gates (docs/numerics.md sections 4 and 7.1): `mpmath` at high precision, every function evaluated at the EXACT doubles it is given.

`mpmath`'s own `gammainc` and `betainc` fail to converge or take minutes for parameters of 1e4 and above, so for those two this file evaluates the defining series and
continued fraction in `mpmath` arithmetic at 100 digits (the working precision grows with the magnitude of the exponents, which cancel). That is an *implementation* of the
same mathematics in different arithmetic, so it is itself gated: `python3 scripts/special_ref.py --selfcheck` compares it with `mpmath`'s native functions on every point where the native
one finishes (relative disagreement must be below 1e-30) and with `scipy.special` at double precision (G3).

Public: ref(fn, args) -> mpf, one entry per function the driver `tools/special-check` implements.
"""
import math
import sys

from mpmath import mp, mpf

mp.dps = 100

SQRT2 = None


def _dps_for(*vals):
    # the exponent of the prefactor is of the size of the parameters; its digits must come from the working precision
    m = max([1.0] + [abs(float(v)) for v in vals if v == v and abs(float(v)) < 1e300])
    return 100 + int(math.log10(m)) + 2


def _lg(x):
    return mp.loggamma(x)


def _gamma_series(a, x, logpre):
    s = mpf(1)
    t = mpf(1)
    n = 0
    eps = mpf(10) ** (-(mp.dps + 5))
    while True:
        n += 1
        t = t * x / (a + n)
        s += t
        if t < eps * s:
            break
        if n > 5_000_000:
            raise RuntimeError("gamma series did not converge")
    return mp.exp(logpre) * s


def _gamma_cf(a, x, logpre):
    # Q = e^-x x^a / Gamma(a) * 1/(x+1-a- 1(1-a)/(x+3-a- 2(2-a)/(x+5-a- ...)))  (modified Lentz)
    tiny = mpf(10) ** (-(mp.dps * 2))
    b = x + 1 - a
    c = 1 / tiny
    d = 1 / b
    h = d
    i = 0
    eps = mpf(10) ** (-(mp.dps + 5))
    while True:
        i += 1
        an = -i * (i - a)
        b += 2
        d = an * d + b
        if d == 0:
            d = tiny
        c = b + an / c
        if c == 0:
            c = tiny
        d = 1 / d
        dl = d * c
        h *= dl
        if abs(dl - 1) < eps:
            break
        if i > 5_000_000:
            raise RuntimeError("gamma continued fraction did not converge")
    return mp.exp(logpre) * a * h


def gamma_pq(a, x):
    """(P, Q) of the regularised incomplete gamma function at exactly a, x (doubles or mpf), both directly accurate."""
    a = mpf(a)
    x = mpf(x)
    if x == 0:
        return mpf(0), mpf(1)
    if x < 1:
        # the series for P converges at once for x < 1; Q = 1 - P is taken at a precision that keeps Q's digits (Q is about a |ln x| for small a)
        extra = 0 if a >= mpf(10) ** -3 else int(-mp.log10(a)) + 10
        with mp.workdps(_dps_for(a, x) + extra):
            a = mpf(a)
            x = mpf(x)
            logpre = a * mp.log(x) - x - _lg(a + 1)
            p = _gamma_series(a, x, logpre)
            q = 1 - p
        return +p, +q
    with mp.workdps(_dps_for(a, x)):
        a = mpf(a)
        x = mpf(x)
        logpre = a * mp.log(x) - x - _lg(a + 1)  # ln(x^a e^-x / Gamma(a+1))
        small = mpf(10) ** -20
        if x < a + 1:
            p = _gamma_series(a, x, logpre)
            q = 1 - p
            if q < small:
                q = _gamma_cf(a, x, logpre)
        else:
            q = _gamma_cf(a, x, logpre)
            p = 1 - q
            if p < small:
                p = _gamma_series(a, x, logpre)
        return +p, +q


def beta_ix(a, b, x):
    """(I_x(a, b), 1 - I_x(a, b)) at exactly a, b, x (0 < x < 1), both directly accurate.

    The cheaper of the two series (by the turning point of the term ratio and its decay rate) gives one tail directly; the other is 1 minus it. When that other is so small that
    the subtraction leaves fewer than 40 correct digits, the same cheap series is rerun at a larger working precision (never the expensive one: for x near 1 and a of 1e7 it
    does not converge in 5e6 terms)."""
    a = mpf(a)
    b = mpf(b)
    x = mpf(x)
    if x <= 0:
        return mpf(0), mpf(1)
    if x >= 1:
        return mpf(1), mpf(0)
    # 1 - x must keep its digits when x is tiny (a double down to 5e-324)
    boost = int(-mp.log10(x)) + 10 if x < mpf(10) ** -3 else 0
    base = _dps_for(a, b) + boost
    for extra in (0, 150, 400, 1200, 3500):
        with mp.workdps(base + extra):
            aa = mpf(a)
            bb = mpf(b)
            xx = mpf(x)
            yy = 1 - xx

            def side(a1, b1, x1, y1):
                # front * 2F1(a1+b1, 1; a1+1; x1) / a1, with x1 below the switch point of the fraction when it is the cheap side
                logfront = a1 * mp.log(x1) + b1 * mp.log(y1) - (_lg(a1) + _lg(b1) - _lg(a1 + b1)) - mp.log(a1)
                sm = mpf(1)
                t = mpf(1)
                n = 0
                eps = mpf(10) ** (-(mp.dps + 5))
                while True:
                    t = t * x1 * (a1 + b1 + n) / (a1 + 1 + n)
                    n += 1
                    sm += t
                    if t < eps * sm:
                        break
                    if n > 3_000_000:
                        raise RuntimeError("beta series did not converge")
                return mp.exp(logfront) * sm

            def cost(a1, b1, x1):
                turn = max(mpf(0), (x1 * (a1 + b1) - a1 - 1) / (1 - x1))
                return turn + (mp.dps * 2.4) / (-mp.log(x1))

            if cost(aa, bb, xx) <= cost(bb, aa, yy):
                lo = side(aa, bb, xx, yy)
                hi = 1 - lo
                direct, other = lo, hi
            else:
                hi = side(bb, aa, yy, xx)
                lo = 1 - hi
                direct, other = hi, lo
            if other == 0 or other > mpf(10) ** -(mp.dps - 45):
                return +lo, +hi
    raise RuntimeError("beta complement below the precision tried")


def lgamma(x):
    with mp.workdps(_dps_for(x) + 20):
        return +mp.loggamma(mpf(x))


def lgamma1p(a):
    """ln Gamma(1 + a) to full relative accuracy for tiny a (the working precision grows with -log10 a)."""
    a = mpf(a)
    with mp.workdps(120 + (int(-mp.log10(a)) if a > 0 and a < 1 else 0)):
        return +mp.loggamma(1 + mpf(a))


def erf(x):
    return +mp.erf(mpf(x))


def erfc(x):
    return +mp.erfc(mpf(x))


def norm_cdf(z):
    return +(mp.erfc(-mpf(z) / mp.sqrt(2)) / 2)


def norm_sf(z):
    return +(mp.erfc(mpf(z) / mp.sqrt(2)) / 2)


def qnorm(p):
    """the z with Phi(z) = p, p an exact double in (0, 1); p > 1/2 by symmetry (1 - p is exact)."""
    p = mpf(p)
    if p == mpf(1) / 2:
        return mpf(0)
    if p > mpf(1) / 2:
        return -qnorm_lower(1 - p)
    return qnorm_lower(p)


def qnorm_lower(p):
    # Newton on f(z) = ln Phi(z) - ln p, started from scipy's double answer, at 100 digits
    from scipy.special import ndtri

    z = mpf(float(ndtri(float(p)))) if float(p) > 0 else mpf(-38.5)
    lp = mp.log(p)
    for _ in range(60):
        cdf = mp.erfc(-z / mp.sqrt(2)) / 2
        f = mp.log(cdf) - lp
        dens = mp.exp(-z * z / 2) / mp.sqrt(2 * mp.pi)
        step = f * cdf / dens
        z = z - step
        if abs(step) < mpf(10) ** (-90) * max(1, abs(z)):
            break
    return +z


def t_tails(t, nu):
    """(cdf, sf) of Student's t with nu degrees of freedom at exactly t, nu."""
    t = mpf(t)
    nu = mpf(nu)
    with mp.workdps(_dps_for(nu) + 10):
        t = mpf(t)
        nu = mpf(nu)
        if t == 0:
            return mpf(1) / 2, mpf(1) / 2
        at = abs(t)
        x = nu / (nu + at * at)
        # two tails of |T| > at is I_x(nu/2, 1/2); one tail is half of it
        two_lo, two_hi = beta_ix(nu / 2, mpf(1) / 2, x)
        tail = two_lo / 2
        other = 1 - tail
        if t > 0:
            return +other, +tail
        return +tail, +other


def chi2_tails(x, nu):
    p, q = gamma_pq(mpf(nu) / 2, mpf(x) / 2)
    return p, q


def f_tails(f, d1, d2):
    """(cdf, sf) of F(d1, d2) at exactly f."""
    f = mpf(f)
    d1 = mpf(d1)
    d2 = mpf(d2)
    if f <= 0:
        return mpf(0), mpf(1)
    with mp.workdps(_dps_for(d1, d2) + 10):
        f = mpf(f)
        d1 = mpf(d1)
        d2 = mpf(d2)
        x = d1 * f / (d1 * f + d2)
        lo, hi = beta_ix(d1 / 2, d2 / 2, x)
        return +lo, +hi


def binom_tails(k, n, p):
    """(P(X <= k), P(X > k)) for X ~ Binomial(n, p), k integer-valued double, 0 <= k < n."""
    k = mpf(k)
    n = mpf(n)
    p = mpf(p)
    # P(X <= k) = I_{1-p}(n - k, k + 1)
    lo, hi = beta_ix(n - k, k + 1, 1 - p)
    return +lo, +hi


# ------------------------------------------------------------------------------------------------------------------------------ quantiles (docs/numerics.md section 9)
# The reference solves tail(x) = tau at the working precision from a START that scipy gives in double precision (only the start: the answer is the root of the 100-digit forward functions above,
# by a secant iteration in ln s that stops at 1e-45 relative). A root below the smallest subnormal is returned as a true tiny number (or 0 when the test below says so), so that the gate
# can apply the same absolute allowance as everywhere else.

TINY = mpf("5e-324")


def _secant_log(g, u0, maxit=80):
    """the root u of the monotone g(u) by the secant method from u0 and u0 + 1e-9"""
    u0 = mpf(u0)
    u1 = u0 + mpf("1e-9")
    f0 = g(u0)
    f1 = g(u1)
    for _ in range(maxit):
        if f1 == f0:
            break
        u2 = u1 - f1 * (u1 - u0) / (f1 - f0)
        if abs(u2 - u1) < mpf("1e-45") * max(mpf(1), abs(u2)):
            return u2
        u0, f0 = u1, f1
        u1 = u2
        f1 = g(u1)
    raise RuntimeError("secant did not converge")


def _bisect_log(g, u0, width=2000, steps=300):
    """the fallback: a sign change by stepping out from u0, then Illinois-style bisection (slow, used only when the start was bad)"""
    u0 = mpf(u0)
    lo, hi = u0 - 1, u0 + 1
    for i in range(24):
        if g(lo) <= 0:
            break
        lo -= 8 * 2 ** i
    for i in range(24):
        if g(hi) >= 0:
            break
        hi += 8 * 2 ** i
    if g(lo) > 0 or g(hi) < 0:
        raise RuntimeError("no sign change")
    for _ in range(steps):
        mid = (lo + hi) / 2
        if g(mid) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < mpf("1e-45") * max(mpf(1), abs(mid)):
            break
    return (lo + hi) / 2


def _solve_log(g, x0):
    """the root, VERIFIED: the residual of the forward function at the answer must be below 1e-35 (the secant method can stop on a flat piece far from the root, which an earlier version of this
    file did once: a quantile whose forward value was 1e-126 for a target of 1e-183 was accepted until the gate's own result, which was right, disagreed)"""
    g0 = g
    lo_u, hi_u = mpf(-746), mpf(740)

    def g(u):
        # every root of these functions is between the smallest subnormal and 1e320, so outside [lo_u, hi_u] the function is continued linearly (steeply, increasing): the iteration cannot wander
        c = min(max(u, lo_u), hi_u)
        return g0(c) + 1000 * (u - c)

    u = None
    try:
        if x0 and x0 > 0 and x0 < 1e308:
            u = _secant_log(g, mp.log(mpf(x0)))
    except RuntimeError:
        u = None
    if u is None or abs(g(u)) > mpf("1e-35"):
        u = _bisect_log(g, mp.log(mpf(x0)) if x0 and x0 > 0 else mpf(0))
    if abs(g(u)) > mpf("1e-35"):
        raise RuntimeError("the root was not found: residual %s" % mp.nstr(g(u), 5))
    return u


def gamma_quantile_ref(tau, a, upper):
    from scipy import special as sp

    tau = mpf(tau)
    a = mpf(a)

    def tail(x):
        p, q = gamma_pq(a, x)
        return q if upper else p

    lt = mp.log(tau)
    t = tail(TINY)
    if (not upper and t >= tau) or (upper and t <= tau):
        return mpf(0)
    sign = 1 if not upper else -1
    g = lambda u: sign * (mp.log(tail(mp.exp(u))) - lt)
    try:
        x0 = float(sp.gammainccinv(float(a), float(tau)) if upper else sp.gammaincinv(float(a), float(tau)))
    except Exception:
        x0 = float(a)
    if not (0.0 < x0 < 1e300):
        x0 = float(a)
    return mp.exp(_solve_log(g, x0))


def beta_pair_ref(tau, a, b, upper):
    """(x, y) with the tail of x equal to tau, the small one solved for, the other 1 minus it (a root below the smallest subnormal answers (0, 1) or (1, 0))"""
    from scipy import special as sp

    tau = mpf(tau)
    a = mpf(a)
    b = mpf(b)
    lt = mp.log(tau)
    # which side is small: the root is at most 1/2 when the lower tail at 1/2 is at least tau (lower target), or the upper tail at 1/2 is at most tau (upper target)
    lo_h, up_h = beta_ix(a, b, mpf(1) / 2)
    x_small = (lo_h >= tau) if not upper else (up_h <= tau)
    # the start in double precision (scipy), used only if it lies on the side found
    try:
        if upper:
            xs = float(sp.betainccinv(float(a), float(b), float(tau)))
        else:
            xs = float(sp.betaincinv(float(a), float(b), float(tau)))
    except Exception:
        xs = float("nan")
    if x_small:
        if not (0.0 < xs <= 0.5):
            xs = 0.25
    else:
        if not (0.5 < xs < 1.0):
            xs = 0.75

    def tails_x(x):
        return beta_ix(a, b, x)

    def tails_y(y):
        lo_y, up_y = beta_ix(b, a, y)
        return up_y, lo_y  # the lower tail of x is the upper tail of y

    def g_of(tails, sign):
        def g(u):
            lo, up = tails(mp.exp(u))
            t = up if upper else lo
            return sign * (mp.log(t) - lt)
        return g

    if x_small:
        lo, up = tails_x(TINY)
        t = up if upper else lo
        if (not upper and t >= tau) or (upper and t <= tau):
            return mpf(0), mpf(1)
        s = mp.exp(_solve_log(g_of(tails_x, 1 if not upper else -1), xs))
        return s, 1 - s
    lo, up = tails_y(TINY)
    t = up if upper else lo
    # for the y side the tail of x DEcreases as y grows
    if (not upper and t <= tau) or (upper and t >= tau):
        return mpf(1), mpf(0)
    ys = 1.0 - xs
    if ys <= 0.0:
        ys = 1e-16
    s = mp.exp(_solve_log(g_of(tails_y, -1 if not upper else 1), ys))
    return 1 - s, s


def t_quantile_ref(p, nu):
    p = mpf(p)
    nu = mpf(nu)
    if p == mpf(1) / 2:
        return mpf(0)
    with mp.workdps(_dps_for(nu) + 10):
        tau = p if p < mpf(1) / 2 else 1 - p
        x, y = beta_pair_ref(2 * tau, nu / 2, mpf(1) / 2, False)
        t = mp.sqrt(nu * y / x) if x > 0 else mp.inf
        return +(t if p > mpf(1) / 2 else -t)


def f_quantile_ref(p, d1, d2, upper):
    d1 = mpf(d1)
    d2 = mpf(d2)
    with mp.workdps(_dps_for(d1, d2) + 10):
        x, y = beta_pair_ref(p, d1 / 2, d2 / 2, upper)
        return +(d2 * x / (d1 * y)) if y > 0 else mp.inf


# name -> (arity, function returning mpf)
def ref(fn, args):
    a = args
    if fn == "lgamma":
        return lgamma(a[0])
    if fn == "lgamma1p":
        return lgamma1p(a[0])
    if fn == "erf":
        return erf(a[0])
    if fn == "erfc":
        return erfc(a[0])
    if fn == "gamma_p":
        return gamma_pq(a[0], a[1])[0]
    if fn == "gamma_q":
        return gamma_pq(a[0], a[1])[1]
    if fn == "beta_lo":
        return beta_ix(a[0], a[1], a[2])[0]
    if fn == "beta_up":
        return beta_ix(a[0], a[1], a[2])[1]
    if fn == "norm_cdf":
        return norm_cdf(a[0])
    if fn == "norm_sf":
        return norm_sf(a[0])
    if fn == "qnorm":
        return qnorm(a[0])
    if fn == "t_cdf":
        return t_tails(a[0], a[1])[0]
    if fn == "t_sf":
        return t_tails(a[0], a[1])[1]
    if fn == "chi2_cdf":
        return chi2_tails(a[0], a[1])[0]
    if fn == "chi2_sf":
        return chi2_tails(a[0], a[1])[1]
    if fn == "f_cdf":
        return f_tails(a[0], a[1], a[2])[0]
    if fn == "f_sf":
        return f_tails(a[0], a[1], a[2])[1]
    if fn == "binom_cdf":
        return binom_tails(a[0], a[1], a[2])[0]
    if fn == "binom_sf":
        return binom_tails(a[0], a[1], a[2])[1]
    if fn == "gamma_quantile":
        return gamma_quantile_ref(a[0], a[1], False)
    if fn == "gamma_quantile_upper":
        return gamma_quantile_ref(a[0], a[1], True)
    if fn == "chi2_quantile":
        return 2 * gamma_quantile_ref(a[0], mpf(a[1]) / 2, False)
    if fn == "chi2_quantile_upper":
        return 2 * gamma_quantile_ref(a[0], mpf(a[1]) / 2, True)
    if fn == "beta_quantile":
        return beta_pair_ref(a[0], a[1], a[2], False)[0]
    if fn == "beta_quantile_comp":
        return beta_pair_ref(a[0], a[1], a[2], False)[1]
    if fn == "beta_quantile_upper":
        return beta_pair_ref(a[0], a[1], a[2], True)[0]
    if fn == "t_quantile":
        return t_quantile_ref(a[0], a[1])
    if fn == "t_quantile_upper":
        return -t_quantile_ref(a[0], a[1])
    if fn == "f_quantile":
        return f_quantile_ref(a[0], a[1], a[2], False)
    if fn == "f_quantile_upper":
        return f_quantile_ref(a[0], a[1], a[2], True)
    raise KeyError(fn)


ARITY = {"lgamma1p": 1, "lgamma": 1, "erf": 1, "erfc": 1, "gamma_p": 2, "gamma_q": 2, "beta_lo": 3, "beta_up": 3, "norm_cdf": 1, "norm_sf": 1, "qnorm": 1,
         "t_cdf": 2, "t_sf": 2, "chi2_cdf": 2, "chi2_sf": 2, "f_cdf": 3, "f_sf": 3, "binom_cdf": 3, "binom_sf": 3,
         "gamma_quantile": 2, "gamma_quantile_upper": 2, "chi2_quantile": 2, "chi2_quantile_upper": 2, "beta_quantile": 3, "beta_quantile_comp": 3, "beta_quantile_upper": 3,
         "t_quantile": 2, "t_quantile_upper": 2, "f_quantile": 3, "f_quantile_upper": 3}


def selfcheck(n=300, seed=1):
    """The reference against mpmath's native functions where they finish and the native value is at least 1e-20 (below that mpmath computes a tail as 1 - the other tail at its own
    precision and returns noise, even a negative number: measured at a = 3.7, b = 578.5, x = 0.578)."""
    import random, signal

    rnd = random.Random(seed)

    class TO(Exception):
        pass

    def alarm(sig, frm):
        raise TO()

    signal.signal(signal.SIGALRM, alarm)
    worst = {"gamma": mpf(0), "beta": mpf(0)}
    done = {"gamma": 0, "beta": 0}
    for _ in range(n):
        a = 10 ** rnd.uniform(-2, 3.5)
        x = a * 10 ** rnd.uniform(-1.5, 0.7)
        p, q = gamma_pq(a, x)
        try:
            signal.alarm(5)
            np_ = mp.gammainc(mpf(a), 0, mpf(x), regularized=True)
            nq = mp.gammainc(mpf(a), mpf(x), mp.inf, regularized=True)
            signal.alarm(0)
        except Exception:
            signal.alarm(0)
            continue
        with mp.workdps(60):
            for r, nat in ((p, np_), (q, nq)):
                if nat > mpf(10) ** -20:
                    worst["gamma"] = max(worst["gamma"], abs(r - nat) / abs(nat))
        done["gamma"] += 1
    for _ in range(n):
        a = 10 ** rnd.uniform(-2, 3)
        b = 10 ** rnd.uniform(-2, 3)
        x = rnd.random()
        lo, hi = beta_ix(a, b, x)
        try:
            signal.alarm(5)
            nl = mp.betainc(mpf(a), mpf(b), 0, mpf(x), regularized=True)
            nh = mp.betainc(mpf(a), mpf(b), mpf(x), 1, regularized=True)
            signal.alarm(0)
        except Exception:
            signal.alarm(0)
            continue
        with mp.workdps(60):
            for r, nat in ((lo, nl), (hi, nh)):
                if nat > mpf(10) ** -20:
                    worst["beta"] = max(worst["beta"], abs(r - nat) / abs(nat))
        done["beta"] += 1
    ok = True
    for k in worst:
        print("selfcheck %s: %d points where mpmath finished, worst relative disagreement %s" % (k, done[k], mp.nstr(worst[k], 3)))
        if done[k] < n // 3 or worst[k] > mpf(10) ** -30:
            ok = False
    return ok


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        sys.exit(0 if selfcheck() else 1)
    print(__doc__)
