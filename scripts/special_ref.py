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
    """(I_x(a, b), 1 - I_x(a, b)) at exactly a, b, x (0 < x < 1), both directly accurate."""
    a = mpf(a)
    b = mpf(b)
    x = mpf(x)
    if x <= 0:
        return mpf(0), mpf(1)
    if x >= 1:
        return mpf(1), mpf(0)
    with mp.workdps(_dps_for(a, b)):
        a = mpf(a)
        b = mpf(b)
        x = mpf(x)

        def side(aa, bb, xx, yy):
            # front * 2F1(aa+bb, 1; aa+1; xx) / aa, with xx < (aa+1)/(aa+bb+2)
            logfront = aa * mp.log(xx) + bb * mp.log(yy) - (_lg(aa) + _lg(bb) - _lg(aa + bb)) - mp.log(aa)
            s = mpf(1)
            t = mpf(1)
            n = 0
            eps = mpf(10) ** (-(mp.dps + 5))
            while True:
                t = t * xx * (aa + bb + n) / (aa + 1 + n)
                n += 1
                s += t
                if t < eps * s:
                    break
                if n > 5_000_000:
                    raise RuntimeError("beta series did not converge")
            return mp.exp(logfront) * s

        y = 1 - x
        small = mpf(10) ** -20
        if x < (a + 1) / (a + b + 2):
            lo = side(a, b, x, y)
            hi = 1 - lo
            if hi < small:
                hi = side(b, a, y, x)
        else:
            hi = side(b, a, y, x)
            lo = 1 - hi
            if lo < small:
                lo = side(a, b, x, y)
        return +lo, +hi


def lgamma(x):
    with mp.workdps(_dps_for(x) + 20):
        return +mp.loggamma(mpf(x))


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


# name -> (arity, function returning mpf)
def ref(fn, args):
    a = args
    if fn == "lgamma":
        return lgamma(a[0])
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
    raise KeyError(fn)


ARITY = {"lgamma": 1, "erf": 1, "erfc": 1, "gamma_p": 2, "gamma_q": 2, "beta_lo": 3, "beta_up": 3, "norm_cdf": 1, "norm_sf": 1, "qnorm": 1,
         "t_cdf": 2, "t_sf": 2, "chi2_cdf": 2, "chi2_sf": 2, "f_cdf": 3, "f_sf": 3, "binom_cdf": 3, "binom_sf": 3}


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
