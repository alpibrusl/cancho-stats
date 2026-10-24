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

    return mpmath.findroot(f, float(ss.t.isf(float(alpha_half), float(df))))


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


def endpoint_ok(got, want, estimate, half, q=0):
    """the interval bound of T-C: absolute 1e-12 * max(|estimate|, half-width), and never below 2 + q units of 2^-1074 (q the quantile of the answer; 2 for t1, which has always been judged by two units): the standard error is a
    double, so in the subnormal range it is off by half a unit, which q multiplies (t2 and paired, whose first run had a case at 2.4 units with q = 5.8)"""
    if got is None:
        return False
    scale = max(abs(mpmath.mpf(estimate)), abs(half))
    floor = (2 + mpmath.mpf(q)) * mpmath.mpf(2) ** -1074   # an endpoint in the subnormal range cannot be closer to the exact value than its spacing
    return abs(mpmath.mpf(got) - want) <= max(mpmath.mpf("1e-12") * scale, floor)


def frac_mpf(fr):
    return mpmath.mpf(fr.numerator) / fr.denominator


def sums(vals):
    f = [Fraction(v) for v in vals]
    return len(f), sum(f, Fraction(0)), sum((v * v for v in f), Fraction(0))


def t2(a_vals, b_vals, pooled, level):
    """two-sample t (Welch, or pooled with `pooled`) of the present values of the groups a and b (docs/test.md section 2), written from the textbook formulas on Fractions, not from the program's grouping of terms:
    m = S1 / n, v = (S2 - S1^2 / n) / (n - 1), t = (m1 - m2) / sqrt(v1/n1 + v2/n2) (Welch) or / sqrt(vp (1/n1 + 1/n2)) (pooled), the Welch-Satterthwaite or n1 + n2 - 2 degrees of freedom, Cohen's d on the pooled sd.
    A dict like t1: n1, n2, statistic, df (Fraction), p (mpf), d, estimate, se, low/high/half (mpf), reason, variance_ratio (a double or None)."""
    n1, s1a, s2a = sums(a_vals)
    n2, s1b, s2b = sums(b_vals)
    out = {"n1": n1, "n2": n2, "statistic": None, "df": None, "p": None, "d": None, "estimate": None, "low": None, "high": None, "half": None, "reason": None, "variance_ratio": None, "se": None, "t2": None}
    if n1 < 2 or n2 < 2:
        out["reason"] = "too-few-rows"
        return out
    m1, m2 = s1a / n1, s1b / n2
    v1 = (s2a - s1a * s1a / n1) / (n1 - 1)
    v2 = (s2b - s1b * s1b / n2) / (n2 - 1)
    delta = m1 - m2
    out["estimate"] = rounded(delta)
    if v2 != 0:
        out["variance_ratio"] = rounded(v1 / v2)
    nu = n1 + n2 - 2
    if pooled:
        vp = ((n1 - 1) * v1 + (n2 - 1) * v2) / nu
        se_sq = vp * (Fraction(1, n1) + Fraction(1, n2))
        df = Fraction(nu)
    else:
        a, b = v1 / n1, v2 / n2
        se_sq = a + b
        df = (a + b) ** 2 / (a * a / (n1 - 1) + b * b / (n2 - 1)) if se_sq != 0 else None
    out["se"] = sqrt_round(se_sq)
    out["df"] = df
    if df is None:
        df_iv = Fraction(1)
    else:
        df_iv = df
    ah = alpha_half(level)
    q = t_quantile_upper(ah, frac_mpf(df_iv))
    se = mpmath.sqrt(frac_mpf(se_sq))
    mean = frac_mpf(delta)
    out["q"] = q
    out["half"] = q * se
    out["low"] = mean - q * se
    out["high"] = mean + q * se
    if se_sq == 0:
        out["reason"] = "variance-zero"
        return out
    tsq = delta * delta / se_sq
    try:
        t = sqrt_round(tsq)
    except OverflowError:
        out["reason"] = "statistic-beyond-double"
        return out
    out["statistic"] = -t if delta < 0 else t
    vp = ((n1 - 1) * v1 + (n2 - 1) * v2) / nu
    d = sqrt_round(delta * delta / vp)
    out["d"] = -d if delta < 0 else d
    out["x"] = df / (df + tsq)
    out["p"] = p_two_sided_t(frac_mpf(df), out["x"])
    out["t2"] = tsq
    return out


def p_two_sided_t(df, x):  # noqa: F811  (a real df: an mpf or an int; x an exact Fraction in [0, 1])
    """the two-sided p-value of Student's t on df degrees of freedom from x = df / (df + t^2): I_x(df / 2, 1 / 2)"""
    if x == 0:
        return mpmath.mpf(0)
    return mpmath.betainc(mpmath.mpf(df) / 2, mpmath.mpf(1) / 2, 0, mpmath.mpf(x.numerator) / x.denominator, regularized=True)


def paired(pairs, level):
    """the paired t: t1 of the exact differences a - b over the pairs where both are present (mu0 = 0)"""
    diffs = [Fraction(a) - Fraction(b) for a, b in pairs]
    return t1(diffs, 0.0, level)


def chi2_sf(df, x2):
    """P(chi-square with df degrees of freedom > x2) at 60 digits, x2 an exact Fraction"""
    if x2 == 0:
        return mpmath.mpf(1)
    return mpmath.gammainc(mpmath.mpf(df) / 2, frac_mpf(x2) / 2, mpmath.inf, regularized=True)


def _chisq_out(n, df, x2, effect_sq, rows, cols, min_exp, below5):
    out = {"n": n, "df": df, "x2": x2, "statistic": None, "p": None, "effect": None, "min_expected": min_exp, "below5": below5, "rows": rows, "cols": cols, "reason": None}
    if n == 0:
        out["reason"] = "too-few-rows"
        return out
    if df < 1:
        out["reason"] = "too-few-levels"
        return out
    out["statistic"] = rounded(x2)
    out["p"] = chi2_sf(df, x2)
    out["effect"] = mpmath.sqrt(frac_mpf(effect_sq))
    return out


def chisq_independence(table):
    """Pearson's chi-square of a table of counts (a list of rows) from the textbook formula on Fractions: expected counts r_i c_j / N, X2 = sum (O - E)^2 / E (the program groups the terms differently)"""
    r = len(table)
    c = len(table[0]) if r else 0
    n = sum(sum(row) for row in table)
    if n == 0:
        return _chisq_out(0, 0, Fraction(0), Fraction(0), 0, 0, None, 0)
    rs = [sum(row) for row in table]
    cs = [sum(table[i][j] for i in range(r)) for j in range(c)]
    x2 = Fraction(0)
    below5 = 0
    min_exp = None
    for i in range(r):
        for j in range(c):
            e = Fraction(rs[i] * cs[j], n)
            x2 += (table[i][j] - e) ** 2 / e
            if e < 5:
                below5 += 1
            if min_exp is None or e < min_exp:
                min_exp = e
    m = min(r, c)
    out = _chisq_out(n, (r - 1) * (c - 1), x2, x2 / (n * max(m - 1, 1)), r, c, min_exp, below5)
    out["expected"] = [[Fraction(rs[i] * cs[j], n) for j in range(c)] for i in range(r)]
    return out


def chisq_fit(counts, props):
    """goodness of fit of the counts to the proportions (Fractions summing to 1): X2 = sum (O - N p)^2 / (N p)"""
    k = len(counts)
    n = sum(counts)
    if n == 0:
        return _chisq_out(0, 0, Fraction(0), Fraction(0), 0, 0, None, 0)
    x2 = Fraction(0)
    below5 = 0
    min_exp = None
    for o, p in zip(counts, props):
        e = n * p
        x2 += (o - e) ** 2 / e
        if e < 5:
            below5 += 1
        if min_exp is None or e < min_exp:
            min_exp = e
    out = _chisq_out(n, k - 1, x2, x2 / n, k, 1, min_exp, below5)
    out["expected"] = [n * p for p in props]
    return out
