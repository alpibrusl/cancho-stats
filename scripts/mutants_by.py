#!/usr/bin/env python3
"""Mutants of --by in tools/stats/main.cho (docs/describe.md section 10): each is applied to a copy, the driver rebuilt, scripts/gate_describe.py rerun; killed = the gate fails or the driver traps. A mutant that does
not compile is INVALID. Exit 1 if one survives or is invalid.   python3 scripts/mutants_by.py [--check]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
F = "tools/stats/main.cho"
F = "tools/stats/main.cho"
M = [
    ("null-group-inverted", F, "if by_null {\n        return missing;", "if by_null {\n        return !missing;"),
    ("missing-cell-matches-key", F, "if missing {\n        return false;\n    }", ""),
    ("key-comparison-inverted", F, "return bytes.equal(line[from..to], key);", "return !bytes.equal(line[from..to], key);"),
    ("quoted-empty-is-missing", F, "let missing = to == from && cs[3 * by + 2] == 0;", "let missing = to == from;"),
    ("filter-skips-column-zero", F, "if by >= 0 && !in_group(cells, line, by, key, by_null) {", "if by > 0 && !in_group(cells, line, by, key, by_null) {"),
    ("filter-ignored-in-rows", F, "if by >= 0 && !in_group(cells, line, by, key, by_null) {\n        return 4;", "if by >= 0 && in_group(cells, line, by, key, by_null) {\n        return 4;"),
    ("not-in-group-is-a-refusal", F, "} else if st != 4 {", "} else if st != 5 {"),
    ("limit-rows-for-cells", F, "if st == 5 {\n        return rule.refuse_int(heap, e, \"limit.too-many-rows\"", "if st == 6 {\n        return rule.refuse_int(heap, e, \"limit.too-many-rows\""),
    ("levels-filter-ignored", F, "mine = scan.in_group(cs, line, by, key, by_null);", "mine = true;"),
    ("levels-row-not-uncounted", F, "if !mine {\n                            rows = rows - 1;", "if !mine {\n                            rows = rows - 0;"),
    ("group-limit-off-by-one", F, "room = map.size(mr) < max_groups;", "room = map.size(mr) <= max_groups;"),
    ("nulls-not-counted", F, "nulls = nulls + 1;", "nulls = nulls + 0;"),
    ("null-group-always-printed", F, "if clean_null && gr.nulls > 0 {", "if clean_null && gr.nulls >= 0 {"),
    ("null-group-not-null", F, 'path.full(tp), by, "", true, true, w, e);', 'path.full(tp), by, "", false, true, w, e);'),
    ("groups-descending", F, "map.key_at(mr1, en)) > 0 {", "map.key_at(mr1, en)) < 0 {"),
    ("ambiguity-not-detected", F, "} else if matches > 1 {", "} else if matches > 2 {"),
    ("no-match-wrong", F, "if matches == 0 {\n                            e = rule.refuse(heap, e, \"stats.no-such-column\"", "if matches == 1 {\n                            e = rule.refuse(heap, e, \"stats.no-such-column\""),
    ("max-groups-ceiling", F, '"max-groups") > 1000 ||', '"max-groups") > 1001 ||'),
    ("max-groups-floor", F, '"max-groups") < 1 {', '"max-groups") < 0 {'),
    ("max-groups-default", F, "|none|50|the most groups", "|none|49|the most groups"),
    ("group-rows-wrong", F, 'w = json.put_key(heap, w, "rows");\n                                                w = json.put_int(heap, w, info.rows);\n                                            } else {', 'w = json.put_key(heap, w, "rows");\n                                                w = json.put_int(heap, w, info.ncols);\n                                            } else {'),
]
sys.path.insert(0, str(ROOT / "scripts"))
import mutsrc  # noqa: E402
M = mutsrc.relocate(ROOT, M)

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
        print("mutants_by: %d mutants, every pattern applies once" % len(M)); return 0
    if not control():
        print("mutants_by: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_by: %d killed of %d, %d survived and explained, %d unexplained" % (killed, len(M), len(M) - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
