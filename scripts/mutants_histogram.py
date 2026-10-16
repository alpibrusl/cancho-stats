#!/usr/bin/env python3
"""Mutants of lib/acc/histogram.cho and of the lines of tools/stats/main.cho that feed it (docs/describe.md S1-G): each is applied to a copy, the driver rebuilt, scripts/gate_describe.py rerun; killed = the gate fails or the driver traps. A mutant that does
not compile is INVALID. Exit 1 if one survives or is invalid.   python3 scripts/mutants_histogram.py [--check]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
M = [
    ("predicate-sign", "lib/acc/histogram.cho", "keys[mid], i - bins, kmin, 0 - i, kmax)", "keys[mid], bins - i, kmin, 0 - i, kmax)"),
    ("predicate-equality", "lib/acc/histogram.cho", "if sign < 0 {\n                lo = mid + 1;", "if sign <= 0 {\n                lo = mid + 1;"),
    ("edge-coefficient", "lib/acc/histogram.cho", "bins - i, kmin, i, kmax, 0, kmin);", "bins - i, kmin, i + 1, kmax, 0, kmin);"),
    ("edge-unit", "lib/acc/histogram.cho", "emin - 1074, sign < 0);", "emin - 1073, sign < 0);"),
    ("search-middle", "lib/acc/histogram.cho", "let mid = (lo + hi) / 2;", "let mid = (lo + hi + 1) / 2;"),
    ("coefficient-sign", "lib/acc/histogram.cho", "cneg != nn);", "cneg == nn);"),
    ("exponent-alignment", "lib/acc/histogram.cho", "contents(wt), 5, ff - emin);", "contents(wt), 5, ff);"),
    ("constant-column", "tools/stats/main.cho", "if kmin == kmax {\n                        nb = 1;\n                    }", ""),
    ("distinct-comparison", "tools/stats/main.cho", "contents(kw)[t] != contents(kw)[t - 1] {", "contents(kw)[t] == contents(kw)[t - 1] {"),
    ("distinct-first", "tools/stats/main.cho", "if t == 0 || contents(kw)[t]", "if t == 1 || contents(kw)[t]"),
    ("counts-difference", "tools/stats/main.cho", "                        prev = upto;", ""),
    ("last-edge", "tools/stats/main.cho", "while e <= nb {", "while e < nb {"),
    ("default-bins", "tools/stats/main.cho", "|none|10|the number of equal-width", "|none|9|the number of equal-width"),
    ("bins-ceiling", "tools/stats/main.cho", 'table, "bins") > 1000 {', 'table, "bins") > 1001 {'),
]

EXPLAINED = {}


import re


def pattern(o):
    """the text of a mutant as a regex that tolerates the indentation `cancho fmt` gives a continuation line"""
    parts = o.split("\n")
    return r"\n\s*".join(re.escape(p.strip()) if i else re.escape(p) for i, p in enumerate(parts))


def apply(t, o, w):
    return re.sub(pattern(o), lambda m: w, t, count=1)


def occurrences(t, o):
    return len(re.findall(pattern(o), t))


def main():
    for n, fl, o, w in M:
        t = (ROOT / fl).read_text()
        if occurrences(t, o) != 1 or o == w:
            print("mutant %s: pattern occurs %d times" % (n, occurrences(t, o)))
            return 1
    if "--check" in sys.argv:
        print("mutants_histogram: %d mutants, every pattern applies once" % len(M)); return 0
    bad, killed = [], 0
    for n, fl, o, w in M:
        with tempfile.TemporaryDirectory() as d:
            d = pathlib.Path(d)
            for x in ("lib", "tools", "generated"): shutil.copytree(ROOT / x, d / x)
            shutil.copy(ROOT / "cancho.toml", d / "cancho.toml")
            (d / "build").mkdir(); shutil.copytree(ROOT / "build" / "deps", d / "build" / "deps")
            f = d / fl; f.write_text(apply(f.read_text(), o, w))
            r = subprocess.run([CANCHO, "build", "--bin", "stats"], cwd=d, capture_output=True, text=True)
            if r.returncode != 0:
                bad.append(n); print("%-28s INVALID: %s" % (n, r.stderr[-120:].replace("\n", " "))); continue
            try:
                g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe.py"), "--driver", str(d / "build" / "stats"), "--files", "80"], capture_output=True, text=True, timeout=300)
                dead = g.returncode != 0
                if not dead:
                    g2 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe_refusals.py"), "--bin", str(d / "build" / "stats"), "--fuzz", "150"], capture_output=True, text=True, timeout=600)
                    dead = g2.returncode != 0
            except subprocess.TimeoutExpired:
                dead = True
        if dead: killed += 1; print("%-28s killed" % n)
        elif n in EXPLAINED: print("%-28s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-28s SURVIVED" % n)
    print("mutants_histogram: %d killed of %d, %d survived and explained, %d unexplained" % (killed, len(M), len(M) - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
