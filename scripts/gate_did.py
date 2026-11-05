#!/usr/bin/env python3
"""S6.1 gates G-D1 to G-D6 for `stats did` (docs/estimation.md sections 2 and 7):

    python3 scripts/gate_did.py [--bin build/stats] [--cases 300] [--seed 1] [--boot-cases 40]

An independent implementation recomputes every answer: the levels are the correctly rounded quotients, the normal equations are built and SOLVED in exact rational arithmetic (fractions), the effect and the naive figure follow
from the definitions, and the bootstrap replicates are recomputed with the draw rule of docs/estimation.md section 2.2 from scripts/oracle_test.py's SplitMix64 (kind 3), each replicate fitted exactly. The program's numbers
must equal the exact ones to the RELATIVE BOUND fixed in docs/estimation.md section 7 (1e-9, against |value| plus a floor of 1e-9 of the largest level so that a coefficient near zero is not judged by its digits); the maximum
gap seen is printed. Also: the same bytes with the rows in another order (the point estimate; the bootstrap depends on the order of the units within an arm, as `ci boot` depends on the order of the rows), numpy's lstsq agrees
to the same bound, and every refusal of section 2.3 is reached with its tag. CONTROLS must be able to fail: a pooled bootstrap (units drawn from all units rather than within their arm), the stream kind 2 in place of 3 and the
percentile positions floor(R (1 + level) / 2) must each leave the answer on at least a tenth of the bootstrap cases.

Exit 0: every number is the expected one. Exit 1: a difference. Exit 2: the harness failed.
"""
import argparse
import json
import pathlib
import subprocess
import sys
import tempfile
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
import oracle_test as ot  # noqa: E402

BOUND = 1e-9


def run(binary, tmp, text, argv):
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "did"] + list(argv) + ["--root", str(tmp), "f.csv"], capture_output=True, text=True, timeout=600)


def table(units, days):
    head = "arm,pre,post" + (",d0,d1" if days else "")
    rows = []
    for u in units:
        rows.append(",".join("" if v is None else str(v) for v in u[:3] + (u[3:] if days else ())))
    return head + "\n" + "\n".join(rows) + "\n"


def levels(units, days, minpre):
    """(x, y, treated) of the eligible units, as doubles, in file order"""
    out = []
    for u in units:
        used = u[:5] if days else u[:3]
        if any(v is None for v in used):
            continue
        arm, pre, post = u[0], float(u[1]), float(u[2])
        if days:
            d0, d1 = float(u[3]), float(u[4])
            if d0 < minpre:
                continue
            x, y = pre / d0, post / d1
        else:
            x, y = pre, post
        out.append((x, y, arm))
    return out


def exact_fit(rows):
    """the exact least squares solution of y = a + b x + tau t and the figures the program prints, or None if singular / one arm"""
    n = len(rows)
    n1 = sum(1 for r in rows if r[2] == 1)
    if n1 == 0 or n1 == n:
        return None
    X = [[Fraction(1), Fraction(r[0]), Fraction(r[2])] for r in rows]
    Y = [Fraction(r[1]) for r in rows]
    A = [[sum(X[i][p] * X[i][q] for i in range(n)) for q in range(3)] for p in range(3)]
    rhs = [sum(X[i][p] * Y[i] for i in range(n)) for p in range(3)]
    M = [A[p][:] + [rhs[p]] for p in range(3)]
    for c in range(3):
        piv = next((r for r in range(c, 3) if M[r][c] != 0), None)
        if piv is None:
            return None
        M[c], M[piv] = M[piv], M[c]
        for r in range(3):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [M[r][k] - f * M[c][k] for k in range(4)]
    a, b, tau = (M[p][3] / M[p][p] for p in range(3))
    mean_t = sum(Y[i] for i in range(n) if rows[i][2] == 1) / n1
    base = mean_t - tau
    sxt = sum(Fraction(r[0]) for r in rows if r[2] == 1)
    return {"a": a, "b": b, "tau": tau, "effect": (tau / base) if base > 0 else None, "naive": (sum(Y[i] for i in range(n) if rows[i][2] == 1) / sxt - 1) if sxt > 0 else None,
            "counterfactual": base}


def close(got, want, scale):
    if want is None:
        return got is None
    if got is None:
        return False
    w = float(want)
    return abs(got - w) <= BOUND * (abs(w) + 1e-9 * scale)


def gap(got, want, scale):
    w = float(want)
    return abs(got - w) / (abs(w) + 1e-9 * scale)


def make_case(rng, big):
    """units as tuples (arm, pre, post, d0, d1): decimals, ties in the levels, a heavy unit sometimes, negative levels sometimes"""
    n = rng.pick([6, 8, 12, 20, 33, 60]) if not big else 200
    days = rng.next() % 2 == 0
    negative = rng.next() % 4 == 0
    heavy = rng.next() % 4 == 0
    units = []
    for i in range(n):
        d0 = rng.pick([14, 20, 28, 35, 60])
        d1 = rng.pick([20, 28, 35, 90])
        base = rng.pick([1, 2, 5, 10, 40]) * (50 if heavy and i == 0 else 1)
        arm = 1 if rng.next() % 3 else 0
        if i == 0:
            arm = 1
        if i == 1:
            arm = 0
        pre = round(base * d0 * (0.5 + rng.u()) * (1 if not negative else (rng.u() * 2 - 0.7)), 2)
        post = round(base * d1 * (0.6 + rng.u()) * (1 + (0.2 if arm else 0.0)) * (1 if not negative else (rng.u() * 2 - 0.5)), 2)
        if not days:
            pre, post = round(pre / d0, 3), round(post / d1, 3)
        units.append((arm, pre, post, d0, d1))
    if n > 6 and rng.next() % 3 == 0:
        i = 2 + rng.next() % (n - 2)
        field = rng.next() % (5 if days else 3)
        if field != 0:
            u = list(units[i])
            u[field] = None
            units[i] = tuple(u)
    return units, days


def parse(r):
    d = json.loads(r.stdout)
    return d["ok"], d.get("data"), d.get("error")


def replicates(rows, seed, reps, kind=3, pooled=False):
    """the exact (tau, effect) of every replicate, in replicate order: treated units then holdout units, each drawn with replacement from its arm (or from all units if pooled)"""
    tr = [i for i, r in enumerate(rows) if r[2] == 1]
    ho = [i for i, r in enumerate(rows) if r[2] == 0]
    n = len(rows)
    out = []
    for b in range(reps):
        picked = []
        for j in range(n):
            if pooled:
                picked.append(ot.boot_index(seed, kind, b, j, n))
            elif j < len(tr):
                picked.append(tr[ot.boot_index(seed, kind, b, j, len(tr))])
            else:
                picked.append(ho[ot.boot_index(seed, kind, b, j, len(ho))])
        f = exact_fit([rows[i] for i in picked])
        out.append(f)
    return out


def interval(reps_fits, lnum, lscale, key):
    vals = sorted(float(f[key]) for f in reps_fits if f is not None and f[key] is not None)
    m = len(vals)
    if m < 100:
        return None
    lo = m * (lscale - lnum) // (2 * lscale)
    hi = (m * (lscale + lnum) + 2 * lscale - 1) // (2 * lscale) - 1
    return vals[lo], vals[hi], lo, hi, m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--cases", type=int, default=300)
    ap.add_argument("--boot-cases", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = SM64(a.seed)
    bad = []
    worst = {"a": 0.0, "b": 0.0, "tau": 0.0, "effect": 0.0, "naive": 0.0}
    numpy_worst = 0.0
    try:
        import numpy as np
    except ImportError:
        np = None
    done = 0
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        # ---- G-D1/G-D2/G-D3: the point estimate, against the exact solution, in two row orders
        for k in range(a.cases):
            units, days = make_case(rng, big=(k % 50 == 49))
            minpre = rng.pick([14, 14, 20, 30])
            extra = ["--pre-days", "d0", "--post-days", "d1", "--min-pre-days", str(minpre)] if days else []
            r = run(a.bin, tmp, table(units, days), ["--arm", "arm", "--pre", "pre", "--post", "post"] + extra)
            ok, data, err = parse(r)
            rows = levels(units, days, minpre)
            want = exact_fit(rows)
            n1 = sum(1 for x in rows if x[2] == 1)
            n0 = len(rows) - n1
            if want is None or len(rows) < 4 or n1 < 2 or n0 < 2:
                if ok:
                    bad.append(("case %d should be refused" % k, data))
                continue
            if not ok:
                if err["rule"] == "stats.did.singular":
                    continue
                bad.append(("case %d refused" % k, err))
                continue
            scale = max(max(abs(x[0]), abs(x[1])) for x in rows)
            c = data["coefficients"]
            for key, got in (("a", c["a"]), ("b", c["b"]), ("tau", c["tau"])):
                worst[key] = max(worst[key], gap(got, want[key], scale))
                if not close(got, want[key], scale):
                    bad.append(("case %d %s" % (k, key), got, float(want[key])))
            for key in ("effect", "naive"):
                if want[key] is None:
                    if data[key] is not None:
                        bad.append(("case %d %s should be null" % (k, key), data[key]))
                else:
                    worst[key] = max(worst[key], gap(data[key], want[key], scale)) if data[key] is not None else 1.0
                    if not close(data[key], want[key], scale):
                        bad.append(("case %d %s" % (k, key), data[key], float(want[key])))
            if not close(data["counterfactual"], want["counterfactual"], scale):
                bad.append(("case %d counterfactual" % k, data["counterfactual"], float(want["counterfactual"])))
            blank = sum(1 for u in units if any(v is None for v in (u[:5] if days else u[:3])))
            if data["missing"] != blank:
                bad.append(("case %d missing" % k, data["missing"], blank))
            if (data["eligible"], data["n_treated"], data["n_holdout"]) != (len(rows), n1, n0):
                bad.append(("case %d counts" % k, data["eligible"], data["n_treated"], data["n_holdout"]))
            if np is not None:
                X = np.array([[1.0, x[0], float(x[2])] for x in rows])
                y = np.array([x[1] for x in rows])
                tau = float(np.linalg.lstsq(X, y, rcond=None)[0][2])
                numpy_worst = max(numpy_worst, abs(tau - c["tau"]) / (abs(tau) + 1e-9 * scale))
            # the rows in the reverse order give the same bytes
            r2 = run(a.bin, tmp, table(list(reversed(units)), days), ["--arm", "arm", "--pre", "pre", "--post", "post"] + extra)
            if r2.stdout != r.stdout:
                bad.append(("case %d order" % k,))
            # the arm coded the other way round and named with --treated-value 0 gives the same bytes
            if k % 4 == 0:
                flipped = [(1 - u[0],) + tuple(u[1:]) for u in units]
                r3 = run(a.bin, tmp, table(flipped, days), ["--arm", "arm", "--pre", "pre", "--post", "post", "--treated-value", "0"] + extra)
                if r3.stdout != r.stdout:
                    bad.append(("case %d treated-value" % k,))
            done += 1
        # ---- G-D4/G-D5: the bootstrap, replicate by replicate
        differs = {"pooled": 0, "kind2": 0, "positions": 0}
        bc = 0
        for k in range(a.boot_cases):
            units, days = make_case(rng, big=False)
            units = [u for u in units]
            extra = ["--pre-days", "d0", "--post-days", "d1", "--min-pre-days", "14"] if days else []
            rows = levels(units, days, 14)
            n1 = sum(1 for x in rows if x[2] == 1)
            if len(rows) < 6 or n1 < 3 or len(rows) - n1 < 3:
                continue
            seed = rng.pick([0, 1, 7, 20240229, 4611686018427387903, 2 ** 40 + 5])
            reps = rng.pick([100, 120, 200, 240])
            lvl = rng.pick(["0.95", "0.9", "0.8", "0.999", "0.5"])
            r = run(a.bin, tmp, table(units, days), ["--arm", "arm", "--pre", "pre", "--post", "post", "--seed", str(seed), "--replicates", str(reps), "--level", lvl] + extra)
            ok, data, err = parse(r)
            if not ok:
                if err["rule"] in ("stats.did.singular", "stats.did.too-few-units"):
                    continue
                bad.append(("boot case %d refused" % k, err))
                continue
            bc += 1
            lnum = int(Fraction(lvl) * 10 ** 3)
            lscale = 10 ** 3
            fits = replicates(rows, seed, reps)
            scale = max(max(abs(x[0]), abs(x[1])) for x in rows)
            iv = data["interval"]
            for key in ("tau", "effect"):
                want = interval(fits, Fraction(lvl).numerator * (lscale // Fraction(lvl).denominator), lscale, key)
                if want is None:
                    if iv.get(key) is not None:
                        bad.append(("boot case %d %s should be null" % (k, key), iv.get(key)))
                    continue
                got = iv[key]
                if got is None:
                    bad.append(("boot case %d %s is null" % (k, key), want))
                    continue
                if (got["positions"]["low"], got["positions"]["high"], got["replicates_defined"]) != (want[2], want[3], want[4]):
                    bad.append(("boot case %d %s positions" % (k, key), got["positions"], want[2:]))
                if not (close(got["low"], want[0], scale) and close(got["high"], want[1], scale)):
                    bad.append(("boot case %d %s endpoints" % (k, key), got["low"], got["high"], want[:2]))
            # controls: the answer must be able to differ from each wrong rule
            for name, wrong in (("pooled", replicates(rows, seed, reps, pooled=True)), ("kind2", replicates(rows, seed, reps, kind=2))):
                w = interval(wrong, Fraction(lvl).numerator * (lscale // Fraction(lvl).denominator), lscale, "tau")
                g = iv.get("tau")
                if w is not None and g is not None and not (close(g["low"], w[0], scale) and close(g["high"], w[1], scale)):
                    differs[name] += 1
            m = sum(1 for f in fits if f is not None)
            wrong_lo, wrong_hi = m * (lscale - lnum) // (2 * lscale), m * (lscale + lnum) // (2 * lscale)
            vals = sorted(float(f["tau"]) for f in fits if f is not None)
            if iv.get("tau") is not None and m >= 100 and (wrong_hi >= m or abs(vals[wrong_hi] - iv["tau"]["high"]) > 0):
                differs["positions"] += 1
        for name, count in differs.items():
            if bc and count * 10 < bc:
                bad.append(("control %s left the answer on %d of %d cases" % (name, count, bc),))
        # the ceiling of --replicates is allowed: 100000 replicates of a table of six units
        r = run(a.bin, tmp, "arm,pre,post\n1,10,16\n1,12,19\n0,11,13\n0,9,12\n0,10,11\n1,8,13\n", ["--arm", "arm", "--pre", "pre", "--post", "post", "--replicates", "100000", "--seed", "1"])
        ok, data, err = parse(r)
        if not ok or data["interval"]["replicates"] != 100000:
            bad.append(("replicates 100000 refused", err))
        # ---- G-D6: every refusal reached
        good = "arm,pre,post\n1,10,16\n1,12,19\n0,11,13\n0,9,12\n0,10,11\n1,8,13\n"
        refusals = [
            ("stats.test-columns", good, ["--arm", "arm", "--pre", "pre"]),
            ("stats.seed-required", good, ["--arm", "arm", "--pre", "pre", "--post", "post", "--replicates", "100"]),
            ("args.bad-value", good, ["--arm", "arm", "--pre", "pre", "--post", "post", "--seed", "x"]),
            ("args.bad-value", good, ["--arm", "arm", "--pre", "pre", "--post", "post", "--replicates", "99", "--seed", "1"]),
            ("args.bad-value", good, ["--arm", "arm", "--pre", "pre", "--post", "post", "--replicates", "100001", "--seed", "1"]),
            ("args.bad-value", good, ["--arm", "arm", "--pre", "pre", "--post", "post", "--treated-value", "2"]),
            ("stats.no-such-column", good, ["--arm", "arm", "--pre", "pre", "--post", "nope"]),
            ("stats.did.arm-not-binary", "arm,pre,post\n1,10,16\n2,12,19\n0,11,13\n0,9,12\n", ["--arm", "arm", "--pre", "pre", "--post", "post"]),
            ("stats.did.no-treated", "arm,pre,post\n0,10,16\n0,12,19\n0,11,13\n0,9,12\n", ["--arm", "arm", "--pre", "pre", "--post", "post"]),
            ("stats.did.no-holdout", "arm,pre,post\n1,10,16\n1,12,19\n1,11,13\n1,9,12\n", ["--arm", "arm", "--pre", "pre", "--post", "post"]),
            ("stats.did.too-few-units", "arm,pre,post\n1,10,16\n0,12,19\n0,11,13\n", ["--arm", "arm", "--pre", "pre", "--post", "post"]),
            ("stats.did.too-few-units", "arm,pre,post\n1,10,16\n1,12,19\n1,11,13\n0,9,12\n", ["--arm", "arm", "--pre", "pre", "--post", "post"]),
            ("stats.did.too-few-units", "arm,pre,post\n1,10,16\n0,12,19\n0,11,13\n0,9,12\n0,8,10\n", ["--arm", "arm", "--pre", "pre", "--post", "post"]),
            ("stats.did.singular", "arm,pre,post\n1,5,16\n1,5,19\n0,5,13\n0,5,12\n0,5,11\n", ["--arm", "arm", "--pre", "pre", "--post", "post"]),
            ("stats.did.days-not-positive", "arm,pre,post,d0,d1\n1,10,16,0,5\n1,12,19,5,5\n0,11,13,5,5\n0,9,12,5,5\n", ["--arm", "arm", "--pre", "pre", "--post", "post", "--pre-days", "d0", "--post-days", "d1", "--min-pre-days", "0"]),
            ("stats.test-columns", good, ["--arm", "arm", "--pre", "pre", "--post", "post", "--pre-days", "pre"]),
        ]
        for tag, text, argv in refusals:
            r = run(a.bin, tmp, text, argv)
            try:
                ok, data, err = parse(r)
            except ValueError:
                bad.append(("refusal %s: no JSON" % tag, r.stdout[:100], r.stderr[:100]))
                continue
            if ok or err["rule"] != tag:
                bad.append(("refusal %s" % tag, ok, err and err["rule"]))
    print("did: %d point cases, %d bootstrap cases; worst relative gap a %.2e b %.2e tau %.2e effect %.2e naive %.2e; numpy tau %.2e; controls differed on %s of %d" % (
        done, bc, worst["a"], worst["b"], worst["tau"], worst["effect"], worst["naive"], numpy_worst, differs, bc))
    if bad:
        for b in bad[:20]:
            print("FAIL", *b)
        print("%d differences" % len(bad))
        return 1
    print("gate_did: ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
