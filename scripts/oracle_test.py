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


def z_quantile(level_text):
    """q > 0 with P(|Z| < q) = level, at 60 digits (the level a decimal text)"""
    num, k = level_parts(level_text)
    return mpmath.sqrt(2) * mpmath.erfinv(mpmath.mpf(num) / mpmath.mpf(10) ** k)


def p_normal_two_sided(z2):
    """2 P(Z > |z|) at 60 digits from the exact z^2 (a Fraction)"""
    return mpmath.erfc(mpmath.sqrt(frac_mpf(z2)) / mpmath.sqrt(2))


def wilson(x, n, q):
    """Wilson's score interval at 60 digits, clamped to [0, 1]"""
    x, n = mpmath.mpf(x), mpmath.mpf(n)
    c = (x + q * q / 2) / (n + q * q)
    h = q * mpmath.sqrt(x * (n - x) / n + q * q / 4) / (n + q * q)
    return max(c - h, mpmath.mpf(0)), min(c + h, mpmath.mpf(1))


def prop_one(x, n, p0, level):
    """one sample: x successes of n against p0 (a Fraction), z^2 = (x - n p0)^2 / (n p0 (1 - p0)) exactly"""
    out = {"n": n, "x": x, "statistic": None, "p": None, "diff": None, "low": None, "high": None, "estimate": None, "reason": None, "z2": None}
    if n == 0:
        out["reason"] = "too-few-rows"
        return out
    z2 = (x - n * p0) ** 2 / (n * p0 * (1 - p0))
    z = sqrt_round(z2)
    out["statistic"] = -z if x < n * p0 else z
    out["diff"] = rounded(Fraction(x, n) - p0)
    out["estimate"] = rounded(Fraction(x, n))
    out["z2"] = z2
    out["p"] = p_normal_two_sided(z2)
    q = z_quantile(level)
    out["q"] = q
    out["low"], out["high"] = wilson(x, n, q)
    out["n_p0"] = rounded(n * p0)
    out["n_q0"] = rounded(n * (1 - p0))
    return out


def prop_two(x1, n1, x2, n2, level):
    """two samples with the pooled z: z^2 = (p1 - p2)^2 / (p (1 - p) (1/n1 + 1/n2)), p = (x1 + x2) / (n1 + n2)"""
    out = {"statistic": None, "p": None, "diff": None, "low": None, "high": None, "h": None, "reason": None, "z2": None, "pooled": None}
    if n1 == 0 or n2 == 0:
        out["reason"] = "too-few-rows"
        return out
    p1, p2 = Fraction(x1, n1), Fraction(x2, n2)
    out["diff"] = rounded(p1 - p2)
    xx, nn = x1 + x2, n1 + n2
    out["pooled"] = rounded(Fraction(xx, nn))
    q = z_quantile(level)
    out["q"] = q
    l1, u1 = wilson(x1, n1, q)
    l2, u2 = wilson(x2, n2, q)
    d = frac_mpf(p1 - p2)
    mp1, mp2 = frac_mpf(p1), frac_mpf(p2)
    out["low"] = d - mpmath.sqrt((mp1 - l1) ** 2 + (u2 - mp2) ** 2)
    out["high"] = d + mpmath.sqrt((u1 - mp1) ** 2 + (mp2 - l2) ** 2)
    out["h"] = 2 * mpmath.asin(mpmath.sqrt(mp1)) - 2 * mpmath.asin(mpmath.sqrt(mp2))
    if xx == 0 or xx == nn:
        out["reason"] = "variance-zero"
        return out
    pp = Fraction(xx, nn)
    z2 = (p1 - p2) ** 2 / (pp * (1 - pp) * (Fraction(1, n1) + Fraction(1, n2)))
    z = sqrt_round(z2)
    out["statistic"] = -z if p1 < p2 else z
    out["z2"] = z2
    out["p"] = p_normal_two_sided(z2)
    return out


def adjust(ps, method):
    """the adjusted p-values of the entries `ps` (floats, None for a test without one) by bonferroni, holm or bh, from the definitions on Fractions of the printed doubles, each rounded once; ties in the order of the entries.
    The sort is Python's (stable) on (p, position); the running maximum and minimum are taken over the exact Fractions (not over rounded values), which is what the definition says."""
    idx = [i for i, p in enumerate(ps) if p is not None]
    m = len(idx)
    exact = {i: Fraction(ps[i]) for i in idx}
    order = sorted(idx, key=lambda i: (exact[i], i))
    out = [None] * len(ps)
    if method == "bonferroni":
        for i in idx:
            out[i] = rounded(min(Fraction(1), m * exact[i]))
    elif method == "holm":
        run = Fraction(0)
        for r, i in enumerate(order, start=1):
            run = max(run, min(Fraction(1), (m - r + 1) * exact[i]))
            out[i] = rounded(run)
    elif method == "bh":
        run = Fraction(1)
        for r in range(m, 0, -1):
            i = order[r - 1]
            run = min(run, min(Fraction(1), m * exact[i] / r))
            out[i] = rounded(run)
    return out


def p_f(d1, d2, f):
    """the upper tail of the F distribution on (d1, d2) degrees of freedom at the exact rational f >= 0 (d2 may be a Fraction): I_x(d2 / 2, d1 / 2) with x = d2 / (d2 + d1 f), 60 digits"""
    if f == 0:
        return mpmath.mpf(1)
    x = Fraction(d2) / (Fraction(d2) + d1 * f)
    return mpmath.betainc(frac_mpf(Fraction(d2)) / 2, mpmath.mpf(d1) / 2, 0, frac_mpf(x), regularized=True)


def anova(groups, welch):
    """one-way analysis of variance of the lists of present values `groups` (docs/rank.md section 3, S4.1), from the textbook formulas on Fractions: SSW = sum (S2_i - S1_i^2 / n_i), SST = S2 - S1^2 / N, SSB = SST - SSW,
    F = (SSB / (k - 1)) / (SSW / (N - k)); Welch's F* = [sum w_i (m_i - mw)^2 / (k - 1)] / [1 + 2 (k - 2) / (k^2 - 1) T] with w_i = n_i / v_i, mw = sum w_i m_i / W, T = sum (1 - w_i / W)^2 / (n_i - 1) and
    df2 = (k^2 - 1) / (3 T), all EXACT rationals. A group with no value is not a group (k counts the groups with at least one). eta^2 = SSB / SST and omega^2 = (SSB - (k - 1) MSW) / (SST + MSW) always from the sums of squares.
    A dict: ns (n of each group as given, 0 for an empty one), k, N, statistic (the double), df1, df2 (Fraction or int), p (mpf), eta, omega (doubles), groups [(n, mean, sd)] for the groups given, reason."""
    out = {"ns": [len(g) for g in groups], "k": 0, "N": 0, "statistic": None, "df1": None, "df2": None, "p": None, "eta": None, "omega": None, "reason": None, "f": None, "groups": []}
    for g in groups:
        n, s1, s2 = sums(g)
        if n == 0:
            out["groups"].append((0, None, None))
            continue
        m = s1 / n
        sd = sqrt_round((s2 - s1 * s1 / n) / (n - 1)) if n >= 2 else None
        out["groups"].append((n, rounded(m), sd))
    live = [(sums(g), g) for g in groups if len(g) > 0]
    k = len(live)
    nn = sum(s[0][0] for s in live)
    out["k"], out["N"] = k, nn
    if k < 2 or nn - k < 1:
        out["reason"] = "too-few-rows"
        return out
    s1 = sum((s[0][1] for s in live), Fraction(0))
    s2 = sum((s[0][2] for s in live), Fraction(0))
    ssw = sum((s[0][2] - s[0][1] * s[0][1] / s[0][0] for s in live), Fraction(0))
    sst = s2 - s1 * s1 / nn
    ssb = sst - ssw
    out["df1"] = k - 1
    if sst != 0:
        out["eta"] = rounded(ssb / sst)
    msw = ssw / (nn - k)
    if sst + msw != 0:
        out["omega"] = rounded((ssb - (k - 1) * msw) / (sst + msw))
    if welch and any(s[0][0] < 2 for s in live):
        out["reason"] = "too-few-rows"   # Welch's weights need a variance in every group; the effect sizes above stand
        return out
    if welch:
        v = [(s[0][2] - s[0][1] * s[0][1] / s[0][0]) / (s[0][0] - 1) for s in live]
        if any(x == 0 for x in v):
            out["reason"] = "variance-zero"
            return out
        # the span of the scaled weights and distances (docs/rank.md section 3): bit lengths of the integers the program divides, a nonzero value more than 360 bits below the largest is `welch-range`
        unit = Fraction(1, 2 ** 1074)
        bw, bd = [], []
        for (n_i, a1, a2), g in live:
            l_int = (n_i * a2 - a1 * a1) / (unit * unit)
            bw.append((n_i * n_i * (n_i - 1)).bit_length() - int(l_int).bit_length())
            dn = abs((a1 / unit) * nn - (s1 / unit) * n_i)
            if dn != 0:
                bd.append(int(dn).bit_length() - (n_i * nn).bit_length())
        out["spread"] = max(max(bw) - min(bw), (max(bd) - min(bd)) if bd else 0)
        w = [Fraction(s[0][0]) / x for s, x in zip(live, v)]
        ww = sum(w, Fraction(0))
        m = [s[0][1] / s[0][0] for s in live]
        mw = sum((wi * mi for wi, mi in zip(w, m)), Fraction(0)) / ww
        t = sum(((1 - wi / ww) ** 2 / (s[0][0] - 1) for wi, s in zip(w, live)), Fraction(0))
        num = sum((wi * (mi - mw) ** 2 for wi, mi in zip(w, m)), Fraction(0)) / (k - 1)
        f = num / (1 + Fraction(2 * (k - 2), k * k - 1) * t)
        out["df2"] = Fraction(k * k - 1) / (3 * t)
    else:
        if ssw == 0:
            out["reason"] = "variance-zero"
            return out
        f = (ssb / (k - 1)) / (ssw / (nn - k))
        out["df2"] = nn - k
    out["f"] = f
    try:
        out["statistic"] = rounded(f)
    except OverflowError:
        out["reason"] = "statistic-beyond-double"
        return out
    out["p"] = p_f(k - 1, out["df2"], f)
    return out
