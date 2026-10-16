#!/usr/bin/env python3
"""Mutants of lib/acc/quantile.cho and of the lines of tools/stats/main.cho that feed it (docs/describe.md S1-G): each is applied to a copy, the driver rebuilt, scripts/gate_describe.py rerun; killed = the gate fails or the driver traps. A mutant that does
not compile is INVALID. Exit 1 if one survives or is invalid.   python3 scripts/mutants_quantile.py [--check]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
M = [
    ("type1-exact-index", "lib/acc/quantile.cho", "var idx = jn + 1;\n            if g == 0 {\n                idx = jn;", "var idx = jn + 1;\n            if g == 0 {\n                idx = jn + 1;"),
    ("type2-average-weights", "lib/acc/quantile.cho", "at_rank(keys, n, jn + 1), 1, 1, 2)", "at_rank(keys, n, jn + 1), 1, 1, 3)"),
    ("type3-offset-sign", "lib/acc/quantile.cho", "let v2 = 2 * np - tenk;", "let v2 = 2 * np + tenk;"),
    ("type3-even-rule", "lib/acc/quantile.cho", "j % 2 == 0 {\n            idx = j;", "j % 2 == 1 {\n            idx = j;"),
    ("type3-negative-guard", "lib/acc/quantile.cho", "if v2 < 0 {", "if v2 < 0 - 1 {"),
    ("type4-b", "lib/acc/quantile.cho", "if qtype == 4 {\n        b24 = 24;", "if qtype == 4 {\n        b24 = 12;"),
    ("type5-a", "lib/acc/quantile.cho", "a24 = 12;\n        b24 = 12;", "a24 = 11;\n        b24 = 12;"),
    ("type7-b", "lib/acc/quantile.cho", "a24 = 24;\n        b24 = 24;", "a24 = 24;\n        b24 = 0;"),
    ("type8-a", "lib/acc/quantile.cho", "a24 = 8;\n        b24 = 8;", "a24 = 9;\n        b24 = 8;"),
    ("type9-b", "lib/acc/quantile.cho", "a24 = 9;\n        b24 = 9;", "a24 = 9;\n        b24 = 8;"),
    ("nppm-formula", "lib/acc/quantile.cho", "pnum * (24 * (n + 1) - a24 - b24);", "pnum * (24 * (n + 1) - a24);"),
    ("fraction-part", "lib/acc/quantile.cho", "let hn = nn - j * den;", "let hn = nn - j * den + 1;"),
    ("blend-weights-swapped", "lib/acc/quantile.cho", "den - hn, hn, den);", "hn, den - hn, den);"),
    ("rank-clamp-low", "lib/acc/quantile.cho", "if i < 1 {", "if i < 0 {"),
    ("rank-clamp-high", "lib/acc/quantile.cho", "if i > n {", "if i > n + 1 {"),
    ("blend-shift", "lib/acc/quantile.cho", "contents(wt), 5, fl - emin);", "contents(wt), 5, fl);"),
    ("blend-sign", "lib/acc/quantile.cho", "contents(ws), len, nl);", "contents(ws), len, nh);"),
    ("blend-unit", "lib/acc/quantile.cho", "emin - 1074, negative);", "emin - 1073, negative);"),
    ("blend-equal-keys", "lib/acc/quantile.cho", "if klo == khi {\n        return (bits_of_key(klo), 0);\n    }", ""),
    ("key-of-value", "tools/stats/main.cho", "rowkeys[j] = rank.key_of(x);", "rowkeys[j] = rank.key_of(0.0);"),
    ("missing-sentinel", "tools/stats/main.cho", "rowkeys[j] = rank.min_int();", "rowkeys[j] = 0;"),
    ("sort-range", "tools/stats/main.cho", "contents(cw), 0, m);\n                let count", "contents(cw), 0, m - 1);\n                let count"),
    ("gather-sentinel", "tools/stats/main.cho", "if k != rank.min_int() {", "if k != 0 {"),
    ("cells-limit-edge", "tools/stats/main.cho", "} else if (info.rows + 1) * info.ncols > max_cells {", "} else if (info.rows + 1) * info.ncols >= max_cells {"),
    ("probability-bound", "tools/stats/main.cho", "num > quantile.pow10(k) {", "num >= quantile.pow10(k) {"),
    ("probability-decimals", "tools/stats/main.cho", "if k > 9 || digits > 12 {", "if k > 10 || digits > 12 {"),
    ("probability-count", "tools/stats/main.cho", "if n > 99 {", "if n > 100 {"),
    ("default-type", "tools/stats/main.cho", "|none|7|which of the nine", "|none|6|which of the nine"),
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
        print("mutants_quantile: %d mutants, every pattern applies once" % len(M)); return 0
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
    print("mutants_quantile: %d killed of %d, %d survived and explained, %d unexplained" % (killed, len(M), len(M) - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
