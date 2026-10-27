#!/usr/bin/env python3
"""Mutants of `--correct` (docs/test.md T-J): lib/stest/correct.cho (the reading of the p-values, the sort, the three corrections, the family object) and the correction parts of tools/stats/main.cho. Each is applied
to a copy, the program rebuilt, scripts/gate_test5.py and the refusals gate rerun; killed = a gate fails or the program traps or hangs. A mutant that does not compile is INVALID. The gates first run on the UNMUTATED
build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.   python3 scripts/mutants_test5.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
C = "lib/stest/correct.cho"
F = "tools/stats/main.cho"
K = "lib/csvx/cat.cho"
M = [
    ("scaled-mantissa-times-divisor", C, "central.limbs_of(m * mult, contents(nw));", "central.limbs_of(m * div, contents(nw));"),
    ("scaled-divisor", C, "central.limbs_of(div, contents(dw));", "central.limbs_of(mult, contents(dw));"),
    ("scaled-exponent", C, "f - 1074, false);", "f - 1073, false);"),
    ("scaled-cap", C, "if float_of_bits(out) > 1.0 {", "if float_of_bits(out) > 2.0 {"),
    ("scaled-subnormal", C, "        f = ef - 1;\n    }\n    let nl = box_slice(heap, 4, 0);", "        f = ef;\n    }\n    let nl = box_slice(heap, 4, 0);"),
    ("sort-direction", C, "while at_ > 0 && float_of_bits(contents(br)[o[at_ - 1]]) > float_of_bits(contents(br)[k]) {", "while at_ > 0 && float_of_bits(contents(br)[o[at_ - 1]]) < float_of_bits(contents(br)[k]) {"),
    ("sort-ties", C, "while at_ > 0 && float_of_bits(contents(br)[o[at_ - 1]]) > float_of_bits(contents(br)[k]) {", "while at_ > 0 && float_of_bits(contents(br)[o[at_ - 1]]) >= float_of_bits(contents(br)[k]) {"),
    ("holm-remaining-tests", C, "scaled(heap, b[o[r - 1]], count - r + 1, 1)", "scaled(heap, b[o[r - 1]], count - r, 1)"),
    ("holm-running-maximum", C, "if float_of_bits(t) > float_of_bits(run) {\n                            run = t;", "if float_of_bits(t) < float_of_bits(run) {\n                            run = t;"),
    ("holm-start", C, "                    var run = 0;\n                    r = 1;", "                    var run = bits_of(1.0);\n                    r = 1;"),
    ("bh-rank", C, "scaled(heap, b[o[r - 1]], count, r)", "scaled(heap, b[o[r - 1]], count, 1)"),
    ("bh-running-minimum", C, "if float_of_bits(t) < float_of_bits(run) {\n                            run = t;", "if float_of_bits(t) > float_of_bits(run) {\n                            run = t;"),
    ("bh-start", C, "var run = bits_of(1.0);\n                    r = count;", "var run = 0;\n                    r = count;"),
    ("bh-from-the-top", C, "                    r = count;\n                    while r >= 1 {", "                    r = count - 1;\n                    while r >= 1 {"),
    ("bonferroni-count", C, "a[o[r - 1]] = scaled(heap, b[o[r - 1]], count, 1);", "a[o[r - 1]] = scaled(heap, b[o[r - 1]], count - 1, 1);"),
    ("method-holm", C, "let is_holm = bytes_eq(method, \"holm\");", "let is_holm = bytes_eq(method, \"hol\");"),
    ("method-bh", C, "let is_bh = bytes_eq(method, \"bh\");", "let is_bh = bytes_eq(method, \"b\");"),
    ("family-size", C, "w = json.put_int(heap, w, count);", "w = json.put_int(heap, w, n);"),
    ("family-null-entries", C, "if contents(lr2)[k] == 1 {\n                    w = json.put_float(heap, w, float_of_bits(contents(ar)[k]));", "if contents(lr2)[k] == 7 {\n                    w = json.put_float(heap, w, float_of_bits(contents(ar)[k]));"),
    ("pattern-length", C, "i = next_p(text, j);", "i = next_p(text, i + 1);"),
    ("max-tests-edge", F, "if named > cli.nat(args, parsed, table, \"max-tests\") {", "if named >= cli.nat(args, parsed, table, \"max-tests\") {"),
    ("max-tests-pairs", F, "named = bytes.count_byte(cli.text(args, parsed, table, \"pairs\"), ',') + 1;", "named = bytes.count_byte(cli.text(args, parsed, table, \"cols\"), ',') + 1;"),
    ("max-tests-ceiling", F, "cli.nat(args, parsed, table, \"max-tests\") > 32640 ||", "cli.nat(args, parsed, table, \"max-tests\") > 32641 ||"),
    ("max-tests-floor", F, "|| cli.nat(args, parsed, table, \"max-tests\") < 1 {", "|| cli.nat(args, parsed, table, \"max-tests\") < 0 {"),
    ("correct-only-for-test", F, "if cli.has(parsed, table, \"correct\") && !is_test {", "if cli.has(parsed, table, \"correct\") && false {"),
    ("family-written", F, "if is_test && !bytes.equal(cli.text(args, parsed, table, \"correct\"), \"none\") {\n                        // (when an error", "if is_test && bytes.equal(cli.text(args, parsed, table, \"correct\"), \"none\") {\n                        // (when an error"),
    ("adjustment-echo", K, "var w = json.put_key(heap, w0, \"adjustment\");\n    w = json.put_string(heap, w, correction);", "var w = json.put_key(heap, w0, \"adjustment\");\n    w = json.put_string(heap, w, \"none\");"),
    ("correct-none-allowed", F, "if is_test && !bytes.equal(cli.text(args, parsed, table, \"correct\"), \"none\") {\n        var named", "if is_test && bytes.equal(cli.text(args, parsed, table, \"correct\"), \"none\") {\n        var named"),
]

EXPLAINED = {
    "pattern-length": "equivalent: the next \"p_value\": after i + 1 is the same as the next after j, the bytes between them hold no such key (a string is escaped)",
    "sort-ties": "equivalent: the order of equal p-values cannot change an adjusted value (the running maximum of Holm and the running minimum of BH make tied entries equal), as docs/test.md section 12 says",
}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_test5.py"), "--bin", binary, "--files", str(files)], 600),
                      ([sys.executable, str(ROOT / "scripts" / "gate_describe_refusals.py"), "--bin", binary, "--fuzz", "150"], 600)):
        try:
            if subprocess.run(argv, capture_output=True, text=True, timeout=tmo).returncode != 0:
                return False
        except subprocess.TimeoutExpired:
            return False
    return True


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
        print("mutants_test5: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 60):
        print("mutants_test5: the unmutated build fails the gates: refusing to count kills")
        return 1
    print("control (unmutated) passes the gates")
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
                bad.append(n); print("%-32s INVALID: %s" % (n, r.stderr[-160:].replace("\n", " "))); continue
            dead = not gates(str(d / "build" / "stats"), 60)
        if dead: killed += 1; print("%-32s killed" % n)
        elif n in EXPLAINED: print("%-32s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-32s SURVIVED" % n)
    total = len(M) if only is None else len([1 for n, f, o, w in M if n in only])
    print("mutants_test5: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
