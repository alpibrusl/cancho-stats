#!/usr/bin/env python3
"""Mutants of --cols and --format text in tools/stats/main.cho (docs/describe.md section 11): each is applied to a copy, the driver rebuilt, scripts/gate_describe.py rerun; killed = the gate fails or the driver traps. A mutant that does
not compile is INVALID. Exit 1 if one survives or is invalid.   python3 scripts/mutants_cols.py [--check]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
F = "tools/stats/main.cho"
F = "tools/stats/main.cho"
F = "tools/stats/main.cho"
M = [
    ("cols-start-all-selected", F, "while j < ncols {\n        sel[j] = 0;", "while j < ncols {\n        sel[j] = 1;"),
    ("cols-empty-item-accepted", F, "if len(tok) == 0 {\n            return (3, i);", "if len(tok) == 99 {\n            return (3, i);"),
    ("cols-ambiguity-not-detected", F, "if matches > 1 {\n            return (2, i);", "if matches > 2 {\n            return (2, i);"),
    ("cols-no-match-wrong", F, "if matches == 0 {\n            return (1, i);", "if matches == 1 {\n            return (1, i);"),
    ("cols-selects-wrong-column", F, "sel[at] = 1;", "sel[0] = 1;"),
    ("cols-name-compare", F, "if bytes.equal(names[from..vec.get(ends, j)], tok) {", "if !bytes.equal(names[from..vec.get(ends, j)], tok) {"),
    ("cols-report-skips-selected", F, "if sel[j] == 1 {\n            var from = 0;", "if sel[j] == 0 {\n            var from = 0;"),
    ("cols-levels-count-unselected", F, "if bad[j] == 1 && sel[j] == 1 &&", "if bad[j] == 1 &&"),
    ("text-mode-ignored", F, 'let text_mode = bytes.equal(cli.text(args, parsed, table, "format"), "text");', "let text_mode = false;"),
    ("flat-index-step", F, "f[6 * (depth - 1) + 1] = f[6 * (depth - 1) + 1] + 1;", "f[6 * (depth - 1) + 1] = f[6 * (depth - 1) + 1] + 2;"),
    ("flat-path-leading-dot", F, "if k > 0 {\n                o = buffer.push(heap, o, byte_of(46));", "if k >= 0 {\n                o = buffer.push(heap, o, byte_of(46));"),
    ("flat-bracket-open", F, "o = buffer.push(heap, o, byte_of(91));", "o = buffer.push(heap, o, byte_of(40));"),
    ("flat-bracket-close", F, "o = buffer.push(heap, o, byte_of(93));", "o = buffer.push(heap, o, byte_of(41));"),
    ("flat-empty-array-not-shown", F, "if f[6 * depth + 4] == 0 {\n                    if c == 125 {", "if f[6 * depth + 4] == 7 {\n                    if c == 125 {"),
    ("flat-empty-is-object", F, 'out = put_line(heap, out, doc, f, depth, "[]");', 'out = put_line(heap, out, doc, f, depth, "{}");'),
    ("flat-string-escape-ignored", F, "if int_of(doc[i]) == 92 {\n                            i = i + 1;\n                        }", ""),
    ("flat-key-expected-after-comma", F, "f[6 * (depth - 1) + 5] = 1;\n                }\n                i = i + 1;", "f[6 * (depth - 1) + 5] = 0;\n                }\n                i = i + 1;"),
    ("flat-literal-end-bracket", F, "&& int_of(doc[i]) != 93 {", "{"),
    ("flat-literal-end-brace", F, "&& int_of(doc[i]) != 125", ""),
    ("flat-key-start", F, "f[6 * top + 2] = start + 1;", "f[6 * top + 2] = start;"),
    ("flat-key-end", F, "f[6 * top + 3] = i - 1;", "f[6 * top + 3] = i;"),
    ("flat-children-not-counted", F, "f[6 * top + 4] = f[6 * top + 4] + 1;\n                    out = put_line", "f[6 * top + 4] = f[6 * top + 4] + 0;\n                    out = put_line"),
    ("flat-container-not-counted", F, "if depth > 0 {\n                    f[6 * (depth - 1) + 4] = f[6 * (depth - 1) + 4] + 1;", "if depth > 99 {\n                    f[6 * (depth - 1) + 4] = f[6 * (depth - 1) + 4] + 1;"),
    ("flat-array-kind", F, "f[6 * depth] = 0;\n                if c == 91 {\n                    f[6 * depth] = 1;", "f[6 * depth] = 0;\n                if c == 91 {\n                    f[6 * depth] = 0;"),
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


def control():
    """the same two gates on the UNMUTATED build must pass, or every kill below would be a gate that cannot run (a flag it does not know exits non-zero too)"""
    b = str(ROOT / "build" / "stats")
    g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe.py"), "--bin", b, "--files", "80"], capture_output=True, text=True, timeout=600)
    g2 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe_refusals.py"), "--bin", b, "--fuzz", "150"], capture_output=True, text=True, timeout=600)
    return g.returncode == 0 and g2.returncode == 0


def main():
    for n, fl, o, w in M:
        t = (ROOT / fl).read_text()
        if occurrences(t, o) != 1 or o == w:
            print("mutant %s: pattern occurs %d times" % (n, occurrences(t, o)))
            return 1
    if "--check" in sys.argv:
        print("mutants_cols: %d mutants, every pattern applies once" % len(M)); return 0
    if not control():
        print("mutants_cols: the unmutated build fails the gates: refusing to count kills")
        return 1
    print("control (unmutated) passes both gates")
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
                g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe.py"), "--bin", str(d / "build" / "stats"), "--files", "80"], capture_output=True, text=True, timeout=300)
                dead = g.returncode != 0
                if not dead:
                    g2 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe_refusals.py"), "--bin", str(d / "build" / "stats"), "--fuzz", "150"], capture_output=True, text=True, timeout=600)
                    dead = g2.returncode != 0
            except subprocess.TimeoutExpired:
                dead = True
        if dead: killed += 1; print("%-28s killed" % n)
        elif n in EXPLAINED: print("%-28s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-28s SURVIVED" % n)
    print("mutants_cols: %d killed of %d, %d survived and explained, %d unexplained" % (killed, len(M), len(M) - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
