#!/usr/bin/env python3
"""Mutants of `stats test prop` (docs/test.md T-J): lib/stest/prop.cho (the exact statistic), lib/stest/answer4.cho (the intervals and the entries), the level registration of lib/csvx/cross.cho that it adds and the prop part
of tools/stats/main.cho. Each is applied to a copy, the program rebuilt, scripts/gate_test4.py, scripts/gate_test3.py (the shared counting pass) and the refusals gate rerun; killed = a gate fails or the program traps or hangs.
A mutant that does not compile is INVALID. The gates first run on the UNMUTATED build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.
    python3 scripts/mutants_test4.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
P = "lib/stest/prop.cho"
A = "lib/stest/answer4.cho"
X = "lib/csvx/cross.cho"
F = "tools/stats/main.cho"
K = "lib/csvx/cat.cho"
M = [
    # one sample
    ("one-no-trial", P, "    if n == 0 {\n        return Prop", "    if n == 1 {\n        return Prop"),
    ("one-numerator-sign", P, "let sd = x * d - n * p;", "let sd = x * d + n * p;"),
    ("one-denominator-q", P, "central.limbs_of(d - p, contents(w3));", "central.limbs_of(p, contents(w3));"),
    ("one-denominator-n", P, "central.limbs_of(n * p, contents(w2));", "central.limbs_of(p, contents(w2));"),
    ("one-difference-denominator", P, "central.limbs_of(n * d, contents(w2));", "central.limbs_of(n, contents(w2));"),
    ("one-statistic-sign", P, "\n                        if sd < 0 && b == 0 {", "\n                        if sd > 0 && b == 0 {"),
    ("one-difference-sign", P, "central.limbs_of(n * d, contents(w2));\n                        let (c, e) = central.ratio_bits(heap, contents(w1), 3, contents(w2), 3, 0, sd < 0);", "central.limbs_of(n * d, contents(w2));\n                        let (c, e) = central.ratio_bits(heap, contents(w1), 3, contents(w2), 3, 0, sd > 0);"),
    # two samples
    ("two-numerator-sign", P, "let sd = x1 * n2 - x2 * n1;", "let sd = x1 * n2 + x2 * n1;"),
    ("two-no-trial", P, "if n1 == 0 || n2 == 0 {", "if n1 == 0 && n2 == 0 {"),
    ("two-variance-zero-one-end", P, "if xx == 0 || xx == nn {", "if xx == 0 {"),
    ("two-n-factor", P, "central.limbs_of(nn, contents(w3));", "central.limbs_of(1, contents(w3));"),
    ("two-pooled-variance", P, "central.limbs_of(xx * (nn - xx), contents(w3));", "central.limbs_of(xx * nn, contents(w3));"),
    ("two-difference-denominator", P, "central.limbs_of(n1 * n2, contents(w2));\n                        let (c, e) = central.ratio_bits(heap, contents(w1), 3, contents(w2), 3, 0, sd < 0);", "central.limbs_of(n1 + n2, contents(w2));\n                        let (c, e) = central.ratio_bits(heap, contents(w1), 3, contents(w2), 3, 0, sd < 0);"),
    ("two-statistic-sign", P, "\n                            if sd < 0 && b == 0 {", "\n                            if sd > 0 && b == 0 {"),
    ("two-second-proportion", P, "let p2 = ratio(heap, x2, n2);", "let p2 = ratio(heap, x1, n1);"),
    # the intervals and the entries
    ("quantile-tail", A, "float_of(scale - num) / float_of(2 * scale)", "float_of(scale - num) / float_of(scale)"),
    ("wilson-centre", A, "let c = (x + q2 / 2.0) / (n + q2);", "let c = (x + q2) / (n + q2);"),
    ("wilson-half-width", A, "x * (n - x) / n + q2 / 4.0", "x * (n - x) / n + q2 / 2.0"),
    ("wilson-clamp", A, "    if lo < 0.0 {\n        lo = 0.0;", "    if lo < 0.01 {\n        lo = 0.0;"),
    ("p-factor-two", A, "var p = 2.0 * sf;", "var p = 1.0 * sf;"),
    ("p-sign-of-z", A, "    if a < 0.0 {\n        a = 0.0 - a;", "    if a < 1.0 {\n        a = 0.0 - a;"),
    ("newcombe-lower", A, "let lo = dd - sqrt((p1 - l1) * (p1 - l1) + (u2 - p2) * (u2 - p2));", "let lo = dd - sqrt((p1 - l1) * (p1 - l1) + (p2 - l2) * (p2 - l2));"),
    ("newcombe-upper", A, "let hi = dd + sqrt((u1 - p1) * (u1 - p1) + (p2 - l2) * (p2 - l2));", "let hi = dd + sqrt((u1 - p1) * (u1 - p1) + (u2 - p2) * (u2 - p2));"),
    ("cohen-h-sign", A, "- 2.0 * math.atan2(sqrt(p2), sqrt(1.0 - p2));", "+ 2.0 * math.atan2(sqrt(p2), sqrt(1.0 - p2));"),
    ("cohen-h-factor", A, "let h = 2.0 * math.atan2(sqrt(p1), sqrt(1.0 - p1))", "let h = 1.0 * math.atan2(sqrt(p1), sqrt(1.0 - p1))"),
    ("effect-name-difference", A, "w = json.put_string(heap, w, \"difference\");", "w = json.put_string(heap, w, \"diff\");"),
    ("effect-name-h", A, "w = json.put_string(heap, w, \"cohen_h\");", "w = json.put_string(heap, w, \"cohen_d\");"),
    ("interval-name", A, "\"difference of proportions\", \"newcombe\"", "\"difference of proportions\", \"wald\""),
    ("levels-swapped", A, "w = json.put_key(heap, w, \"a\");\n    w = json.put_string(heap, w, level_a);", "w = json.put_key(heap, w, \"a\");\n    w = json.put_string(heap, w, level_b);"),
    ("missing-one", A, "w = json.put_int(heap, w, rows - n);", "w = json.put_int(heap, w, rows);"),
    ("missing-two", A, "w = json.put_int(heap, w, rows - n1 - n2);", "w = json.put_int(heap, w, rows - n1);"),
    ("expected-count", A, "float_of(n) * float_of(p) / float_of(d)", "float_of(n) * float_of(d - p) / float_of(d)"),
    ("pooled-proportion", A, "float_of(x1 + x2) / float_of(n1 + n2)", "float_of(x1) / float_of(n1)"),
    ("variance-zero-reason", A, "reason = \"variance-zero\";", "reason = \"too-few-rows\";"),
    ("no-trial-reason", A, "    if r.status == 1 {\n        reason = \"too-few-rows\";\n    } else if r.zc != 0 {\n        reason = \"statistic-beyond-double\";\n    }\n    w = json.put_key(heap, w, \"effect_size\");\n    if r.status == 0 && r.dc == 0 {", "    if r.status == 0 {\n        reason = \"too-few-rows\";\n    } else if r.zc != 0 {\n        reason = \"statistic-beyond-double\";\n    }\n    w = json.put_key(heap, w, \"effect_size\");\n    if r.status == 0 && r.dc == 0 {"),
    # the counting pass
    ("levels-any-needed", X, "if !present && two && levels_any && !(tb == fb", "if !present && two && !(tb == fb"),
    # the verb
    ("success-required", F, "if !cli.has(parsed, table, \"success\") {", "if false {"),
    ("p0-upper-edge", F, "p0n >= quantile.pow10(p0k)", "p0n > quantile.pow10(p0k)"),
    ("p0-lower-edge", F, "!p0ok || p0n == 0 ||", "!p0ok || p0n == 7 ||"),
    ("p0-with-group", F, "cli.has(parsed, table, \"group\") && cli.has(parsed, table, \"p0\") {", "cli.has(parsed, table, \"group\") && !cli.has(parsed, table, \"p0\") {"),
    ("group-levels-two", K, "if two && dr.cols_ != 2 {", "if two && dr.cols_ != 3 {"),
    ("success-row-default", K, "map.get(dr.cross.la, success, 0 - 1)", "map.get(dr.cross.la, success, 0)"),
    ("second-trials", K, "n2 = n2 + contents(dr.counts)[row * dr.cols_ + 1];", "n2 = n2 + contents(dr.counts)[row * dr.cols_];"),
    ("second-successes", K, "x2 = contents(dr.counts)[sr * dr.cols_ + 1];", "x2 = contents(dr.counts)[sr * dr.cols_];"),
    ("p0-scale", K, "let p0 = p0n * (dd / quantile.pow10(p0k));", "let p0 = p0n * (dd / quantile.pow10(p0k + 1));"),
    ("p0-echo", K, "w = json.put_float(heap, w, float_of(p0) / 1000000000.0);", "w = json.put_float(heap, w, float_of(p0) / 100000000.0);"),
    ("prop-dispatch", F, "if is_test && bytes.equal(cli.operand(args, parsed, 1), \"prop\") {", "if false {"),
    ("prop-kind-accepted", F, "&& !bytes.equal(cli.operand(args, parsed, 1), \"prop\") && !bytes.equal(cli.operand(args, parsed, 1), \"anova\") && !bytes.equal(cli.operand(args, parsed, 1), \"kruskal\") && !bytes.equal(cli.operand(args, parsed, 1), \"mannwhitney\") && !bytes.equal(cli.operand(args, parsed, 1), \"wilcoxon\") && !bytes.equal(cli.operand(args, parsed, 1), \"ks\") && !bytes.equal(cli.operand(args, parsed, 1), \"fisher\") {", "&& false {"),
    ("levels-any-passed", K, "cli.text(args, parsed, table, \"group\"), two, true, int_of", "cli.text(args, parsed, table, \"group\"), two, false, int_of"),
    ("columns-after-first", K, "while ok && i < ncols_named {\n        let name = bytes.field(list, ',', i + 1);\n        match place.open_operand(fs, root, rel, full) {\n            Opened::Failed(reason) => {\n                e = fail.io_error(heap, e, reason, false, rel);\n            }\n            Opened::Ok(opened) => {\n                var file = opened;\n                borrow mut file as &!handle in {\n                    let (x, e2) = cross.cross(heap, handle, rel, name, cli.text(args, parsed, table, \"group\")", "while ok && i < 1 {\n        let name = bytes.field(list, ',', i + 1);\n        match place.open_operand(fs, root, rel, full) {\n            Opened::Failed(reason) => {\n                e = fail.io_error(heap, e, reason, false, rel);\n            }\n            Opened::Ok(opened) => {\n                var file = opened;\n                borrow mut file as &!handle in {\n                    let (x, e2) = cross.cross(heap, handle, rel, name, cli.text(args, parsed, table, \"group\")"),
]

EXPLAINED = {}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_test4.py"), "--bin", binary, "--files", str(files)], 400),
                      ([sys.executable, str(ROOT / "scripts" / "gate_test3.py"), "--bin", binary, "--files", "60"], 400),
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
        print("mutants_test4: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 100):
        print("mutants_test4: the unmutated build fails the gates: refusing to count kills")
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
            dead = not gates(str(d / "build" / "stats"), 100)
        if dead: killed += 1; print("%-32s killed" % n)
        elif n in EXPLAINED: print("%-32s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-32s SURVIVED" % n)
    total = len(M) if only is None else len([1 for n, f, o, w in M if n in only])
    print("mutants_test4: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
