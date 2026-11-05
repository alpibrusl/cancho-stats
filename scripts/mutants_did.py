#!/usr/bin/env python3
"""Mutants of `stats did` (docs/estimation.md section 2): lib/csvx/did.cho (the sums, the system, the effect, the draws, the positions, the eligibility, the refusals) and the did parts of tools/stats/main.cho. Each is one single edit
of the source that must make scripts/gate_did.py fail:
    python3 scripts/mutants_did.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
D = "lib/csvx/did.cho"
M = [
    ("draw-kind", D, "rng.derive(rng.derive(seed, 3), b)", "rng.derive(rng.derive(seed, 2), b)"),
    ("draw-replicate", D, "rng.derive(rng.derive(seed, 3), b)", "rng.derive(rng.derive(seed, 3), b + 1)"),
    ("draw-seed", D, "rng.derive(rng.derive(seed, 3), b)", "rng.derive(rng.derive(seed + 1, 3), b)"),
    ("draw-index", D, "wrapping_mul(j + 1, rng.gold())", "wrapping_mul(j, rng.gold())"),
    ("draw-treated-modulus", D, "contents(lr)[rng.uniform(z, n1)]", "contents(lr)[rng.uniform(z, n1 - 1)]"),
    ("draw-holdout-modulus", D, "contents(hr)[rng.uniform(z, n0)]", "contents(hr)[rng.uniform(z, n0 - 1)]"),
    ("draw-all-holdout", D, "                                                                                if j < n1 {", "                                                                                if j < 0 {"),
    ("draw-split", D, "                                                                                if j < n1 {", "                                                                                if j < n1 + 1 {"),
    ("replicates-count", D, "                                                                        while b < reps {", "                                                                        while b < reps - 1 {"),
    ("tau-low-position", D, "let lo_pos = nt * (scale - lnum) / (2 * scale);", "let lo_pos = nt * (scale + lnum) / (2 * scale);"),
    ("tau-high-position", D, "let hi_pos = (nt * (scale + lnum) + 2 * scale - 1) / (2 * scale) - 1;", "let hi_pos = (nt * (scale + lnum) + 2 * scale - 1) / (2 * scale);"),
    ("effect-low-position", D, "let lo_pos2 = ne * (scale - lnum) / (2 * scale);", "let lo_pos2 = ne * (scale + lnum) / (2 * scale);"),
    ("effect-high-position", D, "let hi_pos2 = (ne * (scale + lnum) + 2 * scale - 1) / (2 * scale) - 1;", "let hi_pos2 = (ne * (scale + lnum) + 2 * scale - 1) / (2 * scale);"),
    ("tau-sort", D, "rank.sort_range(contents(tkw), contents(tmw), contents(cnw), 0, nt);", "rank.sort_range(contents(tkw), contents(tmw), contents(cnw), 0, nt - 1);"),
    ("effect-sort", D, "rank.sort_range(contents(ekw), contents(tmw), contents(cnw), 0, ne);", "rank.sort_range(contents(ekw), contents(tmw), contents(cnw), 0, ne - 1);"),
    ("effect-counterfactual", D, "let base = mean_t - r[2];", "let base = mean_t + r[2];"),
    ("effect-mean", D, "let mean_t = syt / f1;", "let mean_t = syt / fn_;"),
    ("naive-offset", D, "out[5] = syt / sxt - 1.0;", "out[5] = syt / sxt - 0.0;"),
    ("matrix-treated-x", D, "                m[5] = sxt;", "                m[5] = sx;"),
    ("matrix-n", D, "                m[0] = fn_;", "                m[0] = f1;"),
    ("matrix-xx", D, "                m[4] = sxx;", "                m[4] = sxt;"),
    ("rhs-treated", D, "                r[2] = syt;", "                r[2] = sy;"),
    ("rhs-xy", D, "                r[1] = sxy;", "                r[1] = sxx;"),
    ("sum-treated-x", D, "                        psum.add(a[oxt..oxt + s1], 1, xb, wk);", "                        psum.add(a[oxt..oxt + s1], 1, yb, wk);"),
    ("sum-cross", D, "xsum.add(a[oxy..oxy + s2], xb, yb, wk);", "xsum.add(a[oxy..oxy + s2], xb, xb, wk);"),
    ("eligibility-edge", D, "if d0 < minpre {", "if d0 <= minpre {"),
    ("level-pre", D, "x = x / d0;", "x = x / d1;"),
    ("level-post", D, "y = y / d1;", "y = y / d0;"),
    ("treated-value", D, "treated = a == 0.0;", "treated = a == 1.0;"),
    ("arm-binary", D, "if a != 0.0 && a != 1.0 {", "if a != 0.0 {"),
    ("units-holdout", D, "kept < 4 || n1 < 2 || n0 < 2", "kept < 4 || n1 < 2"),
    ("units-treated", D, "kept < 4 || n1 < 2 || n0 < 2", "kept < 4 || n0 < 2"),
    ("units-total", D, "kept < 4 || n1 < 2 || n0 < 2", "kept < 3 || n1 < 2 || n0 < 2"),
    ("days-missing", D, "vec.get(sb, i * nc + jd) == rank.min_int() || vec.get(sb, i * nc + je) == rank.min_int()", "vec.get(sb, i * nc + jd) == rank.min_int() && vec.get(sb, i * nc + je) == rank.min_int()"),
    ("days-positive", D, "if d0 <= 0.0 || d1 <= 0.0 {", "if d0 < 0.0 || d1 <= 0.0 {"),
    ("seed-required", D, '    if !cli.has(parsed, table, "seed") {', '    if cli.has(parsed, table, "seed") {'),
    ("replicates-low", D, "if reps < 100 || reps > 100000 {", "if reps < 101 || reps > 100000 {"),
    ("replicates-high", D, "if reps < 100 || reps > 100000 {", "if reps < 100 || reps >= 100000 {"),
    ("treated-value-range", D, 'cli.nat(args, parsed, table, "treated-value") > 1', 'cli.nat(args, parsed, table, "treated-value") > 2'),
    ("days-together", D, 'cli.has(parsed, table, "pre-days") != cli.has(parsed, table, "post-days")', 'cli.has(parsed, table, "pre-days") && cli.has(parsed, table, "post-days")'),
    ("pivot-singular", D, "if best <= 1e-12 * scale {", "if best <= 0.0 {"),
]
EXPLAINED = {
    "pivot-singular": "a design that is singular in exact arithmetic gives a pivot that is exactly zero in the elimination of integer-valued sums",
}


def gates(binary):
    try:
        return subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_did.py"), "--bin", binary, "--cases", "60", "--boot-cases", "10"], capture_output=True, text=True, timeout=900).returncode == 0
    except subprocess.TimeoutExpired:
        return False


def main():
    for n, f, o, w in M:
        t = (ROOT / f).read_text()
        if t.count(o) != 1 or o == w:
            print("mutant %s: pattern occurs %d times in %s" % (n, t.count(o), f))
            return 1
    only = None
    if "--only" in sys.argv:
        only = set(sys.argv[sys.argv.index("--only") + 1].split(","))
    if "--check" in sys.argv:
        print("mutants_did: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats")):
        print("mutants_did: the unmutated build fails the gate: refusing to count kills")
        return 1
    print("control (unmutated) passes the gate")
    bad, killed = [], 0
    for n, f, o, w in M:
        if only is not None and n not in only:
            continue
        with tempfile.TemporaryDirectory() as d:
            d = pathlib.Path(d)
            for x in ("lib", "tools", "generated"): shutil.copytree(ROOT / x, d / x)
            shutil.copy(ROOT / "cancho.toml", d / "cancho.toml")
            (d / "build").mkdir(); shutil.copytree(ROOT / "build" / "deps", d / "build" / "deps")
            fp = d / f; fp.write_text(fp.read_text().replace(o, w, 1))
            r = subprocess.run([CANCHO, "build", "--bin", "stats"], cwd=d, capture_output=True, text=True)
            if r.returncode != 0:
                bad.append(n); print("%-26s INVALID: %s" % (n, r.stderr[-160:].replace("\n", " "))); continue
            dead = not gates(str(d / "build" / "stats"))
        if dead: killed += 1; print("%-26s killed" % n)
        elif n in EXPLAINED: print("%-26s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-26s SURVIVED" % n)
    total = len(M) if only is None else len([1 for n, f, o, w in M if n in only])
    print("mutants_did: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
