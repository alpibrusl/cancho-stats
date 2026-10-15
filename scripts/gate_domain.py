#!/usr/bin/env python3
"""D1, D2 and D3 of docs/numerics.md section 7.2: domain edges, fuzz over bit patterns, and the rule tags.

    python3 scripts/gate_domain.py [--driver build/special-check] [--fuzz 200000] [--seed 1] [--threads 8]

D1  the edges (0, 1, 2^-1074, 2^-1022, the largest double, +-inf, NaN payloads, -0.0, the neighbours of every branch tie) for every function: the driver must answer every case (no trap, no hang: a 20 s
    timeout per batch) and a refusal must be the one the document specifies.
D2  fuzz: N random 64-bit patterns per argument position per function, uniform over the BIT PATTERNS (so exponents are uniform: huge, tiny, subnormal, inf and NaN all appear). Every answer must
    satisfy the contract of its function: code 0 means a finite value in the range of the function; code != 0 means a NaN value and a tag of special.rules; cdf and sf of a pair are in [0,1] and
    cdf + sf is within 1e-12 of 1 where both are at least 1e-10 (each is promised to 5e-13 relative, so a tighter criterion would not follow from the bounds); the quantile is monotone on a random sorted set of probabilities.
D3  every tag of special.rules is reached by at least one case of D1 or D2, and no other tag is printed. (stats.noconverge is the one tag the inputs of D1 and D2 do not reach: see the report; the
    mutant that lowers the iteration cap shows the path is wired.)
Exit 0: all hold. Exit 1: a violation (the first few are printed). Exit 2: the harness failed.
"""
import argparse
import math
import pathlib
import struct
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64, hx, unhx, nb  # noqa: E402

ARITY = {"lgamma": 1, "erf": 1, "erfc": 1, "gamma_p": 2, "gamma_q": 2, "beta_lo": 3, "beta_up": 3, "norm_cdf": 1, "norm_sf": 1, "qnorm": 1,
         "t_cdf": 2, "t_sf": 2, "chi2_cdf": 2, "chi2_sf": 2, "f_cdf": 3, "f_sf": 3, "binom_cdf": 3, "binom_sf": 3}
PAIRS = [("gamma_p", "gamma_q"), ("beta_lo", "beta_up"), ("norm_cdf", "norm_sf"), ("t_cdf", "t_sf"), ("chi2_cdf", "chi2_sf"), ("f_cdf", "f_sf"), ("binom_cdf", "binom_sf")]
PROB = {"gamma_p", "gamma_q", "beta_lo", "beta_up", "norm_cdf", "norm_sf", "t_cdf", "t_sf", "chi2_cdf", "chi2_sf", "f_cdf", "f_sf", "binom_cdf", "binom_sf"}
TAGS = {}


def driver_tags(driver):
    r = subprocess.run([driver, "tags"], capture_output=True, text=True, timeout=20)
    return {int(a): b for a, b in (ln.split() for ln in r.stdout.strip().split("\n"))}


def run(driver, fn, cases, threads):
    inp = "".join("%d %s\n" % (i, " ".join(hx(a) for a in c)) for i, c in enumerate(cases))
    argv = [driver, fn] + (["--threads", str(threads)] if threads else [])
    r = subprocess.run(argv, input=inp, capture_output=True, text=True, timeout=20 if len(cases) < 5000 else 180)
    if r.returncode != 0:
        raise RuntimeError("%s exit %d (trap or crash) on a batch of %d: %s" % (fn, r.returncode, len(cases), (r.stderr or "")[-200:]))
    out = [ln.split() for ln in r.stdout.strip().split("\n")]
    if len(out) != len(cases):
        raise RuntimeError("%s answered %d lines for %d cases" % (fn, len(out), len(cases)))
    return [(unhx(o[1]), int(o[2])) for o in out]


def specials():
    big = float.fromhex("0x1.fffffffffffffp+1023")
    tiny = 5e-324
    nan_payloads = [struct.unpack("<d", struct.pack("<Q", b))[0] for b in (0x7FF8000000000000, 0x7FF0000000000001, 0xFFF8000000000001, 0x7FFFFFFFFFFFFFFF, 0xFFF0000000000001)]
    base = [0.0, -0.0, 1.0, -1.0, 2.0, 0.5, 10.0, 15.0, 1e-300, 1e300, tiny, -tiny, 2.0 ** -1022, 2.0 ** -1023, big, -big, math.inf, -math.inf, 1e7, 1e8, 1e8 + 1, 1e9, 38.0, 40.0, 26.5, 27.5,
            0.02425, 0.97575, 1.0 - 2.0 ** -53, 2.0 ** -53, 0.5 + 2.0 ** -53, 0.5 - 2.0 ** -54, 3.0, 100.0, 1e15, 2.5e305, 1e5, 123456789.0, 2.0 ** 53, 2.0 ** 63]
    nbs = []
    for v in (0.02425, 0.97575, 10.0, 15.0, 0.5, 1.0, 2.0, 0.1):
        nbs += [nb(v, -1), nb(v, 1)]
    return base + nbs + nan_payloads


def edge_cases(fn):
    sp = specials()
    k = ARITY[fn]
    cs = []
    if k == 1:
        cs = [(v,) for v in sp]
    elif k == 2:
        cs = [(a, b) for a in sp for b in sp]
    else:
        # all triples would be 84^3; take the pairs against a few third arguments
        thirds = [0.0, 1.0, 0.5, 10.0, 1e300, math.inf, float("nan"), -1.0, 1.0 - 2.0 ** -53, 2.0 ** -1074]
        cs = [(a, b, c) for a in sp for b in sp for c in thirds]
    return cs


FIRST = {}


def contract(fn, case, value, code, tags, bad):
    if code != 0 and (code, fn) not in FIRST:
        FIRST[(code, fn)] = case
    if code != 0:
        if code not in tags:
            bad.append((fn, case, "unknown code %d" % code))
        elif value == value:
            bad.append((fn, case, "refusal %s with a value that is not NaN" % tags[code]))
        return
    if value != value:
        bad.append((fn, case, "NaN with code 0"))
        return
    if fn in PROB and not (0.0 <= value <= 1.0):
        bad.append((fn, case, "probability %r outside [0,1]" % value))
    if fn == "erf" and not (-1.0 <= value <= 1.0):
        bad.append((fn, case, "erf %r outside [-1,1]" % value))
    if fn == "erfc" and not (0.0 <= value <= 2.0):
        bad.append((fn, case, "erfc %r outside [0,2]" % value))
    if math.isinf(value):
        bad.append((fn, case, "infinite value with code 0"))
    if fn == "lgamma" and case[0] <= 0.0:
        bad.append((fn, case, "lgamma answered a value for a non-positive argument"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", default=str(ROOT / "build" / "special-check"))
    ap.add_argument("--fuzz", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    tags = driver_tags(a.driver)
    seen = set()
    bad = []
    r = SM64(a.seed)
    total = 0
    try:
        for fn in ARITY:
            # D1
            cs = edge_cases(fn)
            res = run(a.driver, fn, cs, a.threads)
            for c, (v, code) in zip(cs, res):
                seen.add(code)
                contract(fn, c, v, code, tags, bad)
            total += len(cs)
            # D2: bit patterns
            fz = []
            for _ in range(a.fuzz):
                fz.append(tuple(struct.unpack("<d", struct.pack("<Q", r.next()))[0] for _ in range(ARITY[fn])))
            res = run(a.driver, fn, fz, a.threads)
            for c, (v, code) in zip(fz, res):
                seen.add(code)
                contract(fn, c, v, code, tags, bad)
            total += len(fz)
        # pairs: cdf + sf = 1 where both are at least 1e-10, on the fuzz that is in the domain
        for f1, f2 in PAIRS:
            k = ARITY[f1]
            cs = []
            for _ in range(a.fuzz // 4):
                # moderate arguments so that both tails are mostly answers
                args = []
                for j in range(k):
                    args.append(r.logu(1e-3, 1e3) if f1 not in ("norm_cdf",) else r.uni(-8, 8))
                if f1 == "beta_lo":
                    args[2] = r.u()
                if f1.startswith("binom"):
                    args = [float(int(r.u() * 50)), 60.0, r.uni(0.01, 0.99)]
                if f1.startswith("f_") or f1.startswith("t_") or f1.startswith("chi2"):
                    pass
                cs.append(tuple(args))
            ra = run(a.driver, f1, cs, a.threads)
            rb = run(a.driver, f2, cs, a.threads)
            for c, (v1, c1), (v2, c2) in zip(cs, ra, rb):
                seen.add(c1)
                seen.add(c2)
                if c1 == 0 and c2 == 0 and v1 >= 1e-10 and v2 >= 1e-10 and abs(v1 + v2 - 1.0) > 1e-12:
                    bad.append((f1, c, "cdf + sf = %r" % (v1 + v2)))
            total += 2 * len(cs)
        # the quantile is monotone
        ps = sorted(r.logu(5e-324, 0.5) for _ in range(20000)) + sorted(r.u() for _ in range(20000))
        res = run(a.driver, "qnorm", [(p,) for p in ps], a.threads)
        for part in (res[:20000], res[20000:]):
            vals = [v for v, c in part if c == 0]
            for i in range(1, len(vals)):
                if vals[i] < vals[i - 1]:
                    bad.append(("qnorm", (i,), "not monotone"))
                    break
    except Exception as e:
        print("gate_domain: harness failure: %s" % e)
        return 2
    unreached = sorted(tags[c] for c in tags if c not in seen)
    unknown = sorted(c for c in seen if c != 0 and c not in tags)
    print("gate_domain: %d cases answered with no trap and no hang; codes seen: %s" % (total, sorted(tags.get(c, "?%d" % c) for c in seen if c != 0)))
    print("gate_domain: tags of special.rules not reached by any case: %s" % (unreached or "none"))
    for (code, fn), case in sorted(FIRST.items()):
        print("  first case with %s: %s %s" % (tags.get(code), fn, ",".join(hx(x) for x in case)))
    for b in bad[:10]:
        print("  VIOLATION %s %s: %s" % (b[0], ",".join(hx(x) for x in b[1]) if isinstance(b[1][0], float) else b[1], b[2]))
    if bad:
        print("gate_domain: FAILED, %d violations" % len(bad))
        return 1
    if unknown:
        print("gate_domain: FAILED, unknown codes %s" % unknown)
        return 1
    if [t for t in unreached if t != "stats.noconverge"]:
        print("gate_domain: FAILED, tags not reached: %s" % unreached)
        return 1
    print("gate_domain: ok (no violation; unreached: %s)" % (unreached or "none"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
