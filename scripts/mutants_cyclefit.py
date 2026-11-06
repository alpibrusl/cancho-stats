#!/usr/bin/env python3
"""Mutants of `stats cycle fit` (docs/estimation.md section 4): lib/csvx/cyclefit.cho (the design, the QR, the covariance block, the Wald statistic, the shrinkage, the refusals). Each is one single edit
of the source that must make scripts/gate_cyclefit.py fail:
    python3 scripts/mutants_cyclefit.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
D = "lib/csvx/cyclefit.cho"
M = [
    ("rank-threshold", D, "if d <= 1e-6 * big {", "if d <= 1e-12 * big {"),
    ("dof", D, "var dof = n - p;", "var dof = n;"),
    ("wald-cross", D, "br * (m00 * br - m01 * ar)", "br * (m00 * br + m01 * ar)"),
    ("shrink-sign", D, "sh = 1.0 - 2.0 / f;", "sh = 1.0 + 2.0 / f;"),
    ("shrink-clamp", D, "if sh < 0.0 {", "if sh < 1.0 {"),
    ("too-short-edge", D, "if n < 3 * (order + 3) {", "if n < 3 * (order + 3) - 1 {"),
    ("angle-period", D, "6.283185307179586 * t[i] / period", "6.283185307179586 * t[i] / (period + 0.001)"),
    ("sin-cos-swap", D, "contents(aw)[i * p + order + 1] = math.sin(ang) * sw;", "contents(aw)[i * p + order + 1] = math.cos(ang) * sw;"),
    ("weight-root", D, "let sw = sqrt(w[i]);", "let sw = w[i];"),
    ("response-weight", D, "contents(yw)[i] = y[i] * sw;", "contents(yw)[i] = y[i];"),
    ("cov-cross", D, "c01 = c01 + q[(p - 2) * p + kk] * q[(p - 1) * p + kk];", "c01 = c01 - q[(p - 2) * p + kk] * q[(p - 1) * p + kk];"),
    ("cov-start", D, "var kk = p - 2;", "var kk = p - 1;"),
    ("rss-start", D, "                            var j = p;\n", "                            var j = p + 1;\n"),
    ("coef-swap", D, "ar = contents(cw)[p - 2];", "ar = contents(cw)[p - 1];"),
    ("duplicate-day", D, "if i > 0 && contents(tsw)[i] == contents(tsw)[i - 1] {", "if i > 100 && contents(tsw)[i] == contents(tsw)[i - 1] {"),
    ("zero-weight", D, "} else if wt > 0.0 {", "} else if wt >= 0.0 {"),
    ("negative-weight", D, "if wt < 0.0 {", "if wt < 0.0 - 1.0 {"),
    ("shrink-flag", D, 'let js = !bytes.equal(cli.text(args, parsed, table, "shrink"), "none");', 'let js = bytes.equal(cli.text(args, parsed, table, "shrink"), "none");'),
    ("order-limit", D, 'cli.nat(args, parsed, table, "order") > 3', 'cli.nat(args, parsed, table, "order") > 4'),
    ("period-positive", D, "if !(period > 0.0) || period > 1e9 {", "if period < 0.0 || period > 1e9 {"),
    ("integer-day", D, "|| d != float_of(truncate(d)) {", "|| d < 0.0 {"),
    ("negative-day", D, "} else if d < 0.0 || d >= 1099511627776.0", "} else if d < 0.0 - 5.0 || d >= 1099511627776.0"),
    ("series-missing", D, "w = json.put_int(heap, w, gn.nulls);", "w = json.put_int(heap, w, 0);"),
    ("days-count", D, "w = json.put_int(heap, w, kept);", "w = json.put_int(heap, w, kept + 1);"),
]
EXPLAINED = {
}


def gates(binary):
    try:
        return subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_cyclefit.py"), "--bin", binary, "--cases", "60"], capture_output=True, text=True, timeout=900).returncode == 0
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
        print("mutants_cyclefit: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats")):
        print("mutants_cyclefit: the unmutated build fails the gate: refusing to count kills")
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
    print("mutants_cyclefit: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
