"""The oracle of `stats test` (docs/test.md section 5): every statistic from `fractions.Fraction` of the doubles the text reads as, written from the formulas of the design and sharing no code with the program;
p-values and quantiles from `mpmath` at 60 digits on the EXACT rational statistic. Imported by scripts/gate_test.py (and by the gates of the later kinds)."""
import math
import sys
import pathlib
from fractions import Fraction

import mpmath

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import gate_describe as gd  # noqa: E402

mpmath.mp.dps = 60
EPS = 2.0 ** -53


def rounded(fr):
    """the correctly rounded double of the exact rational"""
    return gd.to_float(fr)


def sqrt_round(fr):
    return gd.sqrt_round(fr)


def p_two_sided_t(df, x):
    """the two-sided p-value of Student's t on df degrees of freedom from x = df / (df + t^2), an exact rational in [0, 1]: I_x(df / 2, 1 / 2)"""
    if x == 0:
        return mpmath.mpf(0)
    return mpmath.betainc(mpmath.mpf(df) / 2, mpmath.mpf(1) / 2, 0, mpmath.mpf(x.numerator) / x.denominator, regularized=True)


def t_quantile_upper(alpha_half, df):
    """q > 0 with P(T > q) = alpha_half for Student's t on df degrees of freedom (alpha_half an exact Fraction in (0, 1/2]), 60 digits"""
    import scipy.stats as ss
    target = 2 * mpmath.mpf(alpha_half.numerator) / alpha_half.denominator
    d = mpmath.mpf(df)

    def f(q):
        return mpmath.betainc(d / 2, mpmath.mpf(1) / 2, 0, d / (d + q * q), regularized=True) - target

    return mpmath.findroot(f, float(ss.t.isf(float(alpha_half), df)))


def level_parts(text):
    """("0.95") -> (95, 2): the numerator and the number of decimals"""
    whole, _, frac = text.partition(".")
    return int(whole + frac), len(frac)


def alpha_half(text):
    num, k = level_parts(text)
    return Fraction(10 ** k - num, 2 * 10 ** k)


def t1(vals, mu0, level):
    """one-sample t of the present values against mu0 (a float) at the confidence level given as decimal text. A dict: n, statistic, df, p (mpf), d, estimate, low/high (mpf), q, half, reason."""
    n = len(vals)
    out = {"n": n, "statistic": None, "df": None, "p": None, "d": None, "estimate": None, "low": None, "high": None, "half": None, "reason": None, "x": None, "q": None}
    if n == 0:
        out["reason"] = "too-few-rows"
        return out
    f = [Fraction(v) for v in vals]
    s1 = sum(f, Fraction(0))
    s2 = sum((v * v for v in f), Fraction(0))
    out["estimate"] = rounded(s1 / n)
    if n < 2:
        out["reason"] = "too-few-rows"
        return out
    L = n * s2 - s1 * s1
    D = s1 - n * Fraction(mu0)
    out["df"] = n - 1
    se_sq = Fraction(L, 1) / (n * n * (n - 1))
    out["se"] = sqrt_round(se_sq)
    # the interval: estimate -+ q se, the exact estimate and the exact se
    mean = mpmath.mpf(s1.numerator) / s1.denominator / n
    ah = alpha_half(level)
    q = t_quantile_upper(ah, n - 1)
    se = mpmath.sqrt(mpmath.mpf(se_sq.numerator) / se_sq.denominator)
    out["q"] = q
    out["half"] = q * se
    out["low"] = mean - q * se
    out["high"] = mean + q * se
    if L == 0:
        out["reason"] = "variance-zero"
        return out
    t2 = Fraction(n - 1) * D * D / L
    try:
        t = sqrt_round(t2)
    except OverflowError:
        # |t| above the largest double: no statistic, no effect size and no p-value (the p-value is below 3.6e-309 for every df); the interval stands
        out["reason"] = "statistic-beyond-double"
        return out
    out["statistic"] = -t if D < 0 else t
    d2 = Fraction(n - 1) * D * D / (n * L)
    d = sqrt_round(d2)
    out["d"] = -d if D < 0 else d
    x = Fraction(L) / (L + D * D)
    out["x"] = x
    out["p"] = p_two_sided_t(n - 1, x)
    out["t2"] = t2
    return out


def p_ok(got, want, df):
    """the p-value bound of docs/test.md T-B: relative 1e-11 + (df + 1) 2^-53, from 1e-300 up; below it the answer must be under 1e-290; exactly 1 when want is 1"""
    if got is None:
        return False
    if want < mpmath.mpf("1e-300"):
        return got < 1e-290
    return abs(mpmath.mpf(got) - want) <= (mpmath.mpf("1e-11") + (df + 1) * EPS) * want


def endpoint_ok(got, want, estimate, half):
    """the interval bound of T-C: absolute 1e-12 * max(|estimate|, half-width), and never below two units of 2^-1074"""
    if got is None:
        return False
    scale = max(abs(mpmath.mpf(estimate)), abs(half))
    floor = 2 * mpmath.mpf(2) ** -1074   # an endpoint in the subnormal range cannot be closer to the exact value than its spacing
    return abs(mpmath.mpf(got) - want) <= max(mpmath.mpf("1e-12") * scale, floor)
