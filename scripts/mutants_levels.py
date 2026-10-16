#!/usr/bin/env python3
"""Mutants of the levels of a text column in tools/stats/main.cho (docs/describe.md section 9): each is applied to a copy, the driver rebuilt, scripts/gate_describe.py rerun; killed = the gate fails or the driver traps. A mutant that does
not compile is INVALID. Exit 1 if one survives or is invalid.   python3 scripts/mutants_levels.py [--check]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
F = "tools/stats/main.cho"
M = [
    ("key-column-high-byte", F, "byte_of(j / 256));", "byte_of(j / 255));"),
    ("key-column-low-byte", F, "byte_of(j % 256));", "byte_of(j % 255));"),
    ("value-includes-column-byte", F, "key[2..len(key)]", "key[1..len(key)]"),
    ("quoted-empty-is-missing", F, "!(to == from && cs[3 * j + 2] == 0)", "!(to == from)"),
    ("numeric-columns-counted", F, "if bad[j] == 1 && !(to == from", "if !(to == from"),
    ("cap-off-by-one", F, "room = map.size(mr) < max_levels;", "room = map.size(mr) <= max_levels;"),
    ("cap-not-recorded", F, "contents(pw)[j] = 1;", "contents(pw)[j] = 0;"),
    ("distinct-counts-repeats", F, "if !known {\n borrow mut distinct", "if known {\n borrow mut distinct"),
    ("count-step", F, "map.put(heap, mm, key, count + 1)", "map.put(heap, mm, key, count + 2)"),
    ("order-by-count-ascending", F, "if ca != cb {\n        return ca > cb;", "if ca != cb {\n        return ca < cb;"),
    ("tie-by-bytes-descending", F, "map.key_at(m, b)) < 0;\n}\n\nfn column_of", "map.key_at(m, b)) > 0;\n}\n\nfn column_of"),
    ("rank-after-previous", F, "precedes(m, prev, e);", "precedes(m, e, prev);"),
    ("rank-previous-slot", F, "top[col * k + r - 1]", "top[col * k + r]"),
    ("best-first", F, "best[col] < 0 || precedes(m, e, best[col])", "best[col] < 0 || precedes(m, best[col], e)"),
    ("levels-zero-still-counted", F, "if k > 0 && any_text {", "if any_text {"),
    ("levels-default", F, "|none|10|how many of the most frequent levels", "|none|9|how many of the most frequent levels"),
    ("levels-ceiling", F, '"levels") > 100 {', '"levels") > 101 {'),
    ("max-levels-ceiling", F, '"max-levels") > 1000000 ||', '"max-levels") > 1000001 ||'),
    ("max-levels-floor", F, '"max-levels") < 1 {', '"max-levels") < 0 {'),
    ("capped-flag-false", F, 'w = json.put_key(heap, w, "levels_capped");\n        w = json.put_bool(heap, w, true);', 'w = json.put_key(heap, w, "levels_capped");\n        w = json.put_bool(heap, w, false);'),
    ("header-counted-as-row", F, "seen_header = true;", "seen_header = true;\n rows = rows - 1;"),
    ("second-pass-blank-rule", F, "} else if len(line) > 0 || ncols == 1 {\n                    rows = rows + 1;", "} else if len(line) > 0 {\n                    rows = rows + 1;"),
    ("changed-ignored-rows", F, "} else if changed || rows != expected {", "} else if changed {"),
    ("changed-ignored-ragged", F, "if inside || badq || n != ncols || rows > expected {", "if inside || badq {"),
]

EXPLAINED = {
    "changed-ignored-rows": "stats.file-changed is reached only if the file changes between the two passes, which no gate can do deterministically: the case is argued, not tested",
    "changed-ignored-ragged": "same: a record the first pass accepted is refused by the second only if the file changed in between",
}


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
        print("mutants_levels: %d mutants, every pattern applies once" % len(M)); return 0
    if not control():
        print("mutants_levels: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_levels: %d killed of %d, %d survived and explained, %d unexplained" % (killed, len(M), len(M) - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
