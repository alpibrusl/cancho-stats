"""The oracle of `stats cycle detect` (docs/estimation.md sections 1.2 and 3), written from the document and not from any client's code: a periodic component shared by several series, searched on a grid of periods over a
polynomial level, with a permutation test of the best of them. numpy does the least squares; the shuffles are the ones of the document (SplitMix64, scripts/oracle_test.py's `sm_mix`/`sm_derive`), so a p-value is a function of the
data and the seed and can be compared with the program's exactly.

    detect(series, seed, ...) -> dict of the same members as the program's answer

`series` is a list of (days, values, weights) in the ORDER OF THE SERIES KEYS ascending by bytes (the program's order), each already sorted by day, with weights above zero and finite values.
"""
import math

import numpy as np

import oracle_test as ot

GOLD = ot.GOLD


def design(t, period, order):
    cols = [(t / max(float(t.max()), 1.0)) ** k for k in range(order + 1)]
    if period:
        a = 2.0 * np.pi * t / period
        cols += [np.sin(a), np.cos(a)]
    return np.column_stack(cols)


def whitened_q(t, w, period, order):
    q, _ = np.linalg.qr(design(t, period, order) * np.sqrt(w)[:, None])
    return q


def shuffles(y, seed, s, count, kind=4):
    """count shuffles of the values y of series number s: replicate m has the stream derive(derive(derive(seed, 4), s), m); draw j (from 0) picks r = (z >> 1) mod (i + 1) for i = n - 1 - j, z = mix(stream + (j + 1) gold), and swaps
    positions i and r (Fisher-Yates from the end)"""
    n = len(y)
    out = np.empty((count, n))
    for m in range(count):
        sb = ot.sm_derive(ot.sm_derive(ot.sm_derive(seed, kind), s), m)
        v = list(y)
        for j in range(n - 1):
            i = n - 1 - j
            z = ot.sm_mix((sb + (j + 1) * GOLD) & ot.MASK64)
            r = (z >> 1) % (i + 1)
            v[i], v[r] = v[r], v[i]
        out[m] = v
    return out


def detect(series, seed, order=2, min_period=21.0, min_cycles=2.0, oversample=4, permutations=499, alpha=0.004, edge=0.8, kind=4):
    cleaned = []
    for t, y, w in series:
        if len(t) >= 3 * (order + 3):
            cleaned.append((np.asarray(t, float), np.asarray(y, float), np.asarray(w, float)))
    out = {"found": False, "reason": None, "series_used": len(cleaned)}
    if not cleaned:
        out["reason"] = "too-short"
        return out
    lo = min(float(c[0].min()) for c in cleaned)
    span = max(float(c[0].max()) for c in cleaned) - lo + 1.0
    if span / min_cycles < min_period:
        out["reason"] = "too-short"
        return out
    freqs = np.arange(1.0 / (span / min_cycles), 1.0 / min_period + 1e-12, 1.0 / (oversample * span))
    periods = 1.0 / freqs
    num = np.zeros(len(periods))
    ss0 = 0.0
    pnum = np.zeros((len(periods), permutations))
    pss0 = np.zeros(permutations)
    fits = []
    for s, (t0, y, w) in enumerate(cleaned):
        t = t0 - lo
        sw = np.sqrt(w)
        yw = y * sw
        q0 = whitened_q(t, w, None, order)
        base = yw - q0 @ (q0.T @ yw)
        if float(base @ base) <= 0:
            continue
        perm = shuffles(y, seed, s, permutations, kind) * sw
        base_p = np.sum((perm @ q0) ** 2, axis=1)
        pss0 += np.sum(perm * perm, axis=1) - base_p
        base_fit = float(np.sum((q0.T @ yw) ** 2))
        ss0 += float(base @ base)
        for i, per in enumerate(periods):
            q = whitened_q(t, w, per, order)
            num[i] += float(np.sum((q.T @ yw) ** 2)) - base_fit
            pnum[i] += np.sum((perm @ q) ** 2, axis=1) - base_p
        fits.append((t, w, yw, base_fit))
    if ss0 <= 0:
        out["reason"] = "too-short"
        return out
    gains = num / ss0
    best = int(np.argmax(gains))
    top = (pnum / np.maximum(pss0, 1e-300)).max(axis=0)
    p = float((1 + np.sum(top >= gains[best])) / (permutations + 1))
    out.update(grid=[float(x) for x in periods], gains=[float(x) for x in gains], best_period=float(periods[best]), best_gain=float(gains[best]), p=p)
    if p > alpha:
        out["reason"] = "p-above-alpha"
        return out
    if periods[best] > edge * periods.max():
        out["reason"] = "edge"
        return out

    def gain(per):
        return sum(float(np.sum((whitened_q(t, w, per, order).T @ yw) ** 2)) - bf for t, w, yw, bf in fits) / ss0
    lo_p, hi_p = periods[min(best + 1, len(periods) - 1)], periods[max(best - 1, 0)]
    for _ in range(30):
        m1, m2 = lo_p + (hi_p - lo_p) / 3.0, hi_p - (hi_p - lo_p) / 3.0
        if gain(m1) < gain(m2):
            lo_p = m1
        else:
            hi_p = m2
    period = 0.5 * (lo_p + hi_p)
    amps = []
    for t, w, yw, _ in fits:
        coef = np.linalg.lstsq(design(t, period, order) * np.sqrt(w)[:, None], yw, rcond=None)[0]
        amps.append(float(np.hypot(coef[-2], coef[-1])))
    out.update(found=True, period=float(period), amplitude=float(np.sqrt(np.mean(np.square(amps)))), gain=float(gains[best]))
    return out
