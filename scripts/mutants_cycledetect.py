#!/usr/bin/env python3
"""Mutants of `stats cycle detect` (docs/estimation.md section 3): lib/csvx/cycledetect.cho (the shuffles, the grid, the statistic, the p-value, the decision, the refinement, the flags). Each is one single edit
of the source that must make scripts/gate_cycledetect.py fail:
    python3 scripts/mutants_cyclefit.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
D = "lib/csvx/cycledetect.cho"
M = [
    ("shuffle-kind", D, "rng.derive(rng.derive(rng.derive(seed, 4), s), m)", "rng.derive(rng.derive(rng.derive(seed, 5), s), m)"),
    ("shuffle-series", D, "rng.derive(rng.derive(rng.derive(seed, 4), s), m)", "rng.derive(rng.derive(rng.derive(seed, 4), s + 1), m)"),
    ("shuffle-replicate", D, "rng.derive(rng.derive(rng.derive(seed, 4), s), m)", "rng.derive(rng.derive(rng.derive(seed, 4), s), m + 1)"),
    ("shuffle-bound", D, "let r = rng.uniform(z, i + 1);", "let r = rng.uniform(z, i);"),
    ("shuffle-draws", D, "    while j < n - 1 {\n        let i = n - 1 - j;", "    while j < n - 2 {\n        let i = n - 1 - j;"),
    ("grid-start-count", D, "    let start = 1.0 / (span / minc);\n    let step = 1.0 / (float_of(over) * span);\n    let stop", "    let start = 1.0 / (span / minc) + 1e-4;\n    let step = 1.0 / (float_of(over) * span);\n    let stop"),
    ("grid-start-core", D, "    let start = 1.0 / (span / minc);\n    let step = 1.0 / (float_of(over) * span);\n    let periods", "    let start = 1.0 / (span / minc) + 1e-4;\n    let step = 1.0 / (float_of(over) * span);\n    let periods"),
    ("grid-step", D, "    let start = 1.0 / (span / minc);\n    let step = 1.0 / (float_of(over) * span);\n    let periods", "    let start = 1.0 / (span / minc);\n    let step = 1.0 / (float_of(over) * span) * 1.01;\n    let periods"),
    ("grid-period", D, "contents(prw)[k] = 1.0 / (start + float_of(k) * step);", "contents(prw)[k] = 1.0 / (start + float_of(k + 1) * step);"),
    ("gain-divisor", D, "contents(nur)[k] / ss0;", "contents(nur)[k] / (ss0 + 1.0);"),
    ("p-count", D, "p = float_of(1 + count) / float_of(perms + 1);", "p = float_of(count) / float_of(perms + 1);"),
    ("p-denominator", D, "p = float_of(1 + count) / float_of(perms + 1);", "p = float_of(1 + count) / float_of(perms);"),
    ("shuffle-pss0", D, "pss0[m] = pss0[m] + sq - bp;", "pss0[m] = pss0[m] + sq;"),
    ("period-num", D, "num[pi] = num[pi] + proj2(contents(qw), n, p1, yw, o) - bf;", "num[pi] = num[pi] + proj2(contents(qw), n, p1, yw, o);"),
    ("period-pnum", D, "pnum[pi * perms + m] + proj2(contents(qw), n, p1, contents(pmw), m * n) - contents(bpw)[m];", "pnum[pi * perms + m] + proj2(contents(qw), n, p1, contents(pmw), m * n);"),
    ("shuffle-weight", D, "contents(pmw)[m * n + i] = contents(pmw)[m * n + i] * sqrt(ww[o + i]);", "contents(pmw)[m * n + i] = contents(pmw)[m * n + i];"),
    ("value-weight", D, "yw[o + i] = yy[o + i] * sqrt(ww[o + i]);", "yw[o + i] = yy[o + i];"),
    ("design-angle", D, "let ang = 6.283185307179586 * t[o + i] / per;", "let ang = 6.283185307179586 * t[o + i] / (per + 0.01);"),
    ("ternary-direction", D, "                if g1 < g2 {", "                if g1 > g2 {"),
    ("ternary-start", D, "lo_p = contents(pr3)[ib];", "lo_p = contents(pr3)[ia];"),
    ("ternary-steps", D, "while it < 30 {", "while it < 20 {"),
    ("amplitude", D, "amp = sqrt(sumsq / float_of(nused));", "amp = sumsq / float_of(nused);"),
    ("alpha-test", D, "if p > alpha {", "if p > alpha * 2.0 {"),
    ("edge-test", D, "} else if contents(pr2)[best] > edge * contents(pr2)[0] {", "} else if contents(pr2)[best] > edge * contents(pr2)[1] {"),
    ("need-count", D, "} else if kept >= need {\n                                used = used + 1;", "} else if kept > need {\n                                used = used + 1;"),
    ("need-fill", D, "} else if kept >= need {\n                                                    contents(ofw)[s] = at;", "} else if kept > need {\n                                                    contents(ofw)[s] = at;"),
    ("shift-days", D, "contents(ttw)[at + i] = contents(tr2)[i] - lo;", "contents(ttw)[at + i] = contents(tr2)[i];"),
    ("dropped", D, "w = json.put_int(heap, w, dropped);", "w = json.put_int(heap, w, 0);"),
    ("found-flag", D, "w = json.put_bool(heap, w, code == 0);", "w = json.put_bool(heap, w, code == 1);"),
    ("seed-required", D, '    if !cli.has(parsed, table, "seed") {', '    if cli.has(parsed, table, "seed") {'),
    ("shuffles-low", D, "if perms < 99 || perms > 9999 {", "if perms < 100 || perms > 9999 {"),
    ("alpha-range", D, "if !ok3 || !(alpha > 0.0) || alpha >= 1.0 {", "if !ok3 || !(alpha > 0.0) || alpha > 1.0 {"),
    ("oversample-low", D, "if over < 1 || over > 16 {", "if over < 0 || over > 16 {"),
    ("series-flag", D, 'if !cli.has(parsed, table, "day") || !cli.has(parsed, table, "y") || !cli.has(parsed, table, "series") {', 'if !cli.has(parsed, table, "day") || !cli.has(parsed, table, "y") {'),
]
EXPLAINED = {
    "shift-days": "equivalent: the gain, the amplitude and the p-value do not depend on where the days start (a shift moves the phase of the sine and cosine and the basis of the polynomial level, not the spaces they span); the shift is kept for the conditioning of the level",
}


def gates(binary):
    try:
        return subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_cycledetect.py"), "--bin", binary, "--cases", "40", "--no-controls"], capture_output=True, text=True, timeout=900).returncode == 0
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
        print("mutants_cycledetect: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats")):
        print("mutants_cycledetect: the unmutated build fails the gate: refusing to count kills")
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
    print("mutants_cycledetect: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
