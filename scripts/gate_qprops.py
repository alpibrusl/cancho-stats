#!/usr/bin/env python3
"""The properties of the quantile functions that need no reference (docs/numerics.md 9.4 item 3), run through `tools/special-check`:

    python3 scripts/gate_qprops.py [--driver build/special-check] [--n 3000] [--seed 1]

  symmetry   t_quantile(p) = -t_quantile(1 - p) bit for bit for p in [1/2, 1); t_quantile_upper(q) = -t_quantile(q); t_quantile(1/2) = 0 exactly.
  scale      chi2_quantile(p, nu) = 2 gamma_quantile(p, nu/2) bit for bit (both tails).
  split      beta_quantile + beta_quantile_comp = 1 to 2.3e-16 (x and y = 1 - x, each answered on its own), and the median of the symmetric beta (a = b) is 1/2 within the bound of 9.2 (the reported number is the largest error as a fraction of that bound).
  reciprocal f_quantile(p, d1, d2) f_quantile_upper(p, d2, d1) = 1 within the bound of 9.2 (the lower quantile of F(d1, d2) is the reciprocal of the upper one of F(d2, d1)).
  iterations the largest iteration count of the sweeps' generators must stay at most 80 (the cap is 200): a start or a derivative that is wrong shows here before it shows as a refusal.
  round trip the forward function at the answer equals the probability: |F(x_hat) / p - 1| reported (gamma, beta, t, F; p from 1e-6 to 0.9), and required below 1e-9 (the forward bound times the
             reciprocal elasticity, loose on purpose: the reference gate is the one that decides the accuracy).
Exit 0: every property holds. Exit 1: one fails. Exit 2: the harness failed.
"""
import argparse
import math
import pathlib
import struct
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64, hx, unhx, bits  # noqa: E402


ITER_MARGIN = 80


def run(driver, fn, cases):
    inp = "".join("%d %s\n" % (i, " ".join(hx(a) for a in c)) for i, c in enumerate(cases))
    r = subprocess.run([driver, fn], input=inp, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        raise RuntimeError("%s exit %d" % (fn, r.returncode))
    out = [ln.split() for ln in r.stdout.strip().split("\n")]
    return [(unhx(o[1]), int(o[2])) for o in out]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", default=str(ROOT / "build" / "special-check"))
    ap.add_argument("--n", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    stats = {}
    n = a.n

    def prob_lower():
        return r.logu(1e-12, 0.5)

    try:
        # symmetry of t
        nus = [r.logu(0.1, 1e6) for _ in range(n)]
        ps = [0.5 + r.u() * 0.4999999 for _ in range(n)]
        up = run(a.driver, "t_quantile", [(p, nu) for p, nu in zip(ps, nus)])
        dn = run(a.driver, "t_quantile", [(1.0 - p, nu) for p, nu in zip(ps, nus)])
        for p, nu, (u, cu), (d, cd) in zip(ps, nus, up, dn):
            if cu == 0 and cd == 0 and bits(u) != bits(-d):
                bad.append("t symmetry p=%r nu=%r: %r vs %r" % (p, nu, u, -d))
        qu = run(a.driver, "t_quantile_upper", [(1.0 - p, nu) for p, nu in zip(ps, nus)])
        for p, nu, (u, cu), (v, cv) in zip(ps, nus, up, qu):
            if cu == 0 and cv == 0 and bits(u) != bits(v):
                bad.append("t upper p=%r nu=%r: %r vs %r" % (p, nu, u, v))
        (m, cm), = run(a.driver, "t_quantile", [(0.5, 7.0)])
        if bits(m) != 0 or cm != 0:
            bad.append("t_quantile(1/2) is %r" % m)
        # scale of chi-square
        cs = [(prob_lower() if i & 1 else 1.0 - prob_lower() * 0.5, r.logu(0.02, 1e5)) for i in range(n)]
        for fn_c, fn_g in (("chi2_quantile", "gamma_quantile"), ("chi2_quantile_upper", "gamma_quantile_upper")):
            rc = run(a.driver, fn_c, cs)
            rg = run(a.driver, fn_g, [(p, nu / 2.0) for p, nu in cs])
            for c, (x, cx), (g, cg) in zip(cs, rc, rg):
                if cx == 0 and cg == 0 and bits(x) != bits(2.0 * g):
                    bad.append("%s %r: %r vs 2 * %r" % (fn_c, c, x, g))
        # split of the beta
        bs = [(prob_lower() if i & 1 else r.uni(0.01, 0.99), r.logu(0.05, 500.0), r.logu(0.05, 500.0)) for i in range(n)]
        bx = run(a.driver, "beta_quantile", bs)
        by = run(a.driver, "beta_quantile_comp", bs)
        worst = 0.0
        for c, (x, cx), (y, cy) in zip(bs, bx, by):
            if cx == 0 and cy == 0:
                worst = max(worst, abs(x + y - 1.0))
                if abs(x + y - 1.0) > 2.3e-16:
                    bad.append("beta split %r: x + y = %r" % (c, x + y))
        stats["split"] = worst
        sym = [(0.5, v, v) for v in (r.logu(0.05, 1e4) for _ in range(300))]
        worst = 0.0
        for c, (x, cx) in zip(sym, run(a.driver, "beta_quantile", sym)):
            kappa = 1.0 / min(1.0, c[1])
            tol = 0.5 * kappa * (1e-13 + 1e-15 * math.log(2.0))
            worst = max(worst, abs(x - 0.5) / tol)
            if abs(x - 0.5) > tol:
                bad.append("median of beta(%r, %r) is %r" % (c[1], c[2], x))
        stats["median"] = worst
        # reciprocal of F
        fs = [(r.uni(0.001, 0.999), r.logu(0.5, 1e4), r.logu(0.5, 1e4)) for _ in range(n)]
        lo = run(a.driver, "f_quantile", fs)
        hi = run(a.driver, "f_quantile_upper", [(p, d2, d1) for p, d1, d2 in fs])
        worst = 0.0
        for c, (u, cu), (v, cv) in zip(fs, lo, hi):
            if cu == 0 and cv == 0:
                e = abs(u * v - 1.0)
                worst = max(worst, e)
                kappa = 1.0 / min(1.0, c[1] / 2, c[2] / 2)
                if e > 4.0 * kappa * (1e-13 + 1e-15 * abs(math.log(c[0]))) * 2 + 1e-15:
                    bad.append("F reciprocal %r: %r" % (c, u * v))
        stats["F reciprocal"] = worst
        # round trip
        rt = {}
        specs = [("gamma", "gamma_quantile", "gamma_p", lambda: (r.logu(1e-6, 0.9), r.logu(0.05, 1e4)), lambda p, a_, x: (a_, x)),
                 ("t", "t_quantile", "t_cdf", lambda: (r.uni(1e-6, 0.9), r.logu(0.3, 1e5)), lambda p, a_, x: (x, a_)),
                 ("beta", "beta_quantile", "beta_lo", lambda: (r.logu(1e-6, 0.9), r.logu(0.05, 300.0), r.logu(0.05, 300.0)), lambda p, a_, b, x: (a_, b, x)),
                 ("F", "f_quantile", "f_cdf", lambda: (r.uni(1e-6, 0.9), r.logu(0.5, 1e4), r.logu(0.5, 1e4)), lambda p, a_, b, x: (x, a_, b))]
        for name, qfn, ffn, gen, mk in specs:
            cs = [gen() for _ in range(n // 2)]
            qs = run(a.driver, qfn, cs)
            fw = []
            keep = []
            for c, (x, code) in zip(cs, qs):
                # a beta answer above 0.9 is skipped: the forward function takes x alone, and 1 - x is then not exact (the pair `beta_quantile` and `beta_quantile_comp` is the way to carry it)
                if code == 0 and x > 0 and x <= 0.9 and name == "beta":
                    keep.append(c)
                    fw.append(mk(*c, x))
                elif code == 0 and x > 0 and name != "beta":
                    keep.append(c)
                    fw.append(mk(*c, x))
            res = run(a.driver, ffn, fw)
            worst = 0.0
            for c, (v, code) in zip(keep, res):
                if code == 0 and c[0] > 0:
                    e = abs(v / c[0] - 1.0)
                    worst = max(worst, e)
                    if e > 1e-9:
                        bad.append("round trip %s %r: F(x_hat) / p = %r" % (name, c, v / c[0]))
            rt[name] = worst
        stats["round trip"] = rt
        # the iteration counts: the cap is 200, and the margin is a measurement (docs/numerics.md 9.3)
        import gate_special as gs
        mx = {}
        for name, fn, gen in (("gamma", "gamma_iterations", gs.g_q_gamma), ("gamma big", "gamma_iterations", gs.g_q_gamma_big), ("beta", "beta_iterations", gs.g_q_beta), ("beta big", "beta_iterations", gs.g_q_beta_big)):
            rr = gs.SM64(77)
            cs = [gen(rr, i) for i in range(n // 2)]
            res = run(a.driver, fn, cs)
            worst = max([v for v, code in res if code == 0] + [0])
            mx[name] = int(worst)
            if worst > ITER_MARGIN:
                bad.append("iterations %s: %d exceeds %d" % (name, worst, ITER_MARGIN))
        stats["iterations"] = mx
    except Exception as e:
        print("gate_qprops: harness failure: %s" % e)
        return 2
    print("gate_qprops: %d cases per property; |x + y - 1| max %.3g; median of the symmetric beta: largest error / bound %.3g; F reciprocal max |u v - 1| %.3g; round trip max |F(x_hat)/p - 1|: %s; largest iteration count (cap 200, required at most %d): %s" % (
        n, stats["split"], stats["median"], stats["F reciprocal"], {k: float("%.3g" % v) for k, v in stats["round trip"].items()}, ITER_MARGIN, stats["iterations"]))
    for b in bad[:8]:
        print("  " + b)
    print("gate_qprops: %s" % ("FAILED, %d violations" % len(bad) if bad else "ok"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
