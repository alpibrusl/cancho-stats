#!/usr/bin/env python3
"""Mutants of `stats ci prop` (docs/resample.md section 3): lib/stest/answerci.cho (Clopper-Pearson, the entries), lib/csvx/ci.cho (the counting pass and the method rules), the ci parts of tools/stats/main.cho. Each is applied to a copy, the program rebuilt, scripts/gate_ci_prop.py and the refusals gate rerun; killed = a gate fails or the program traps or hangs. A mutant that does not compile is INVALID. The gates first run on the UNMUTATED build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.
    python3 scripts/mutants_ci.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
W = "lib/stest/answerci.cho"
C = "lib/csvx/ci.cho"
R = "lib/csvx/run.cho"
F = "tools/stats/main.cho"
M = [
    ("half-alpha", W, "return float_of(scale - num) / float_of(2 * scale);", "return float_of(scale - num) / float_of(scale);"),
    ("cp-lower-shape", W, "quantile.beta_quantile(a2, float_of(x), float_of(n - x + 1))", "quantile.beta_quantile(a2, float_of(x), float_of(n - x))"),
    ("cp-upper-shape", W, "quantile.beta_quantile_upper(a2, float_of(x + 1), float_of(n - x))", "quantile.beta_quantile_upper(a2, float_of(x), float_of(n - x))"),
    ("cp-upper-tail", W, "quantile.beta_quantile_upper(a2, float_of(x + 1), float_of(n - x))", "quantile.beta_quantile(a2, float_of(x + 1), float_of(n - x))"),
    ("cp-lower-guard", W, "        if x > 0 {\n            let (q1, k1)", "        if x > 1 {\n            let (q1, k1)"),
    ("cp-upper-guard", W, "        if x < n {\n            let (q2, k2)", "        if x < n - 1 {\n            let (q2, k2)"),
    ("cp-code", W, "if c1 == 0 && c2 == 0 {", "if c1 == 0 || c2 == 0 {"),
    ("cp-method", W, "\"proportion\", \"clopper-pearson\", level", "\"proportion\", \"wilson\", level"),
    ("branch", W, "        if exact {\n            w = answer4.put_interval(heap, w, \"proportion\", \"clopper-pearson\"", "        if !exact {\n            w = answer4.put_interval(heap, w, \"proportion\", \"clopper-pearson\""),
    ("estimate", W, "    if r.status == 0 {\n        w = json.put_float(heap, w, float_of_bits(r.p1b));", "    if r.status == 1 {\n        w = json.put_float(heap, w, float_of_bits(r.p1b));"),
    ("too-few", W, "    if r.status == 1 {\n        reason = \"too-few-rows\";\n    } else if exact {", "    if r.status == 2 {\n        reason = \"too-few-rows\";\n    } else if exact {"),
    ("missing-one", W, "w = json.put_int(heap, w, x);\n    w = json.put_key(heap, w, \"missing\");\n    w = json.put_int(heap, w, rows - n);", "w = json.put_int(heap, w, x);\n    w = json.put_key(heap, w, \"missing\");\n    w = json.put_int(heap, w, rows);"),
    ("missing-two", W, "w = json.put_int(heap, w, rows - n1 - n2);", "w = json.put_int(heap, w, rows - n1);"),
    ("levels", W, "w = json.put_string(heap, w, level_a);", "w = json.put_string(heap, w, level_b);"),
    ("newcombe-low", W, "lo = dd - sqrt((p1 - l1) * (p1 - l1) + (u2 - p2) * (u2 - p2));", "lo = dd - sqrt((p1 - l1) * (p1 - l1) + (p2 - l2) * (p2 - l2));"),
    ("newcombe-high", W, "hi = dd + sqrt((u1 - p1) * (u1 - p1) + (p2 - l2) * (p2 - l2));", "hi = dd + sqrt((u1 - p1) * (u1 - p1) + (u2 - p2) * (u2 - p2));"),
    ("newcombe-method", W, "\"difference of proportions\", \"newcombe\"", "\"difference of proportions\", \"wilson\""),
    ("exact-flag", C, "if bytes.equal(m, \"clopper-pearson\") {\n            exact = true;", "if bytes.equal(m, \"wilson\") {\n            exact = true;"),
    ("two-methods", C, "if two && !bytes.equal(m, \"newcombe\") {", "if two && false {"),
    ("one-methods", C, "} else if !two && !bytes.equal(m, \"wilson\") && !exact {", "} else if !two && false {"),
    ("empty-item", C, "if len(bytes.field(list, ',', i + 1)) == 0 {", "if len(bytes.field(list, ',', i + 1)) == 1 {"),
    ("success-row", C, "let sr = map.get(dr.cross.la, success, 0 - 1);", "let sr = map.get(dr.cross.la, success, 0);"),
    ("n-one", C, "n1 = n1 + contents(dr.counts)[row * dr.cols_];", "n1 = n1 + 1;"),
    ("n-two", C, "n2 = n2 + contents(dr.counts)[row * dr.cols_ + 1];", "n2 = n2 + 0;"),
    ("x-two", C, "x2 = contents(dr.counts)[sr * dr.cols_ + 1];", "x2 = contents(dr.counts)[sr * dr.cols_];"),
    ("group-levels", C, "if two && dr.cols_ != 2 {", "if two && dr.cols_ != 3 {"),
    ("header", C, "w = json.put_string(heap, w, \"prop\");", "w = json.put_string(heap, w, \"mean\");"),
    ("level-echo", C, "w = json.put_key(heap, w, \"level\");\n        w = json.put_float(heap, w, level);", "w = json.put_key(heap, w, \"level\");\n        w = json.put_float(heap, w, 0.95);"),
    ("verb", F, "} else if bytes.equal(cli.operand(args, parsed, 0), \"ci\") {\n        is_ci = true;", "} else if bytes.equal(cli.operand(args, parsed, 0), \"cix\") {\n        is_ci = true;"),
    ("file-operand", F, "is_ci = true;\n        fidx = 2;", "is_ci = true;\n        fidx = 1;"),
    ("kind", F, "} else if !bytes.equal(cli.operand(args, parsed, 1), \"prop\") && !bytes.equal(cli.operand(args, parsed, 1), \"mean\") && !bytes.equal(cli.operand(args, parsed, 1), \"diff\") && !bytes.equal(cli.operand(args, parsed, 1), \"corr\") && !bytes.equal(cli.operand(args, parsed, 1), \"boot\") {", "} else if !bytes.equal(cli.operand(args, parsed, 1), \"prop\") && !bytes.equal(cli.operand(args, parsed, 1), \"mean\") && !bytes.equal(cli.operand(args, parsed, 1), \"diff\") && false && !bytes.equal(cli.operand(args, parsed, 1), \"boot\") {"),
    ("operands", F, "} else if cli.operand_count(parsed) != 3 {", "} else if cli.operand_count(parsed) > 3 {"),
    ("required-flags", F, "bytes.equal(cli.operand(args, parsed, 1), \"prop\") && (!cli.has(parsed, table, \"cols\") || !cli.has(parsed, table, \"success\")) {", "bytes.equal(cli.operand(args, parsed, 1), \"prop\") && (!cli.has(parsed, table, \"cols\") && !cli.has(parsed, table, \"success\")) {"),
    ("dispatch", F, "if is_ci {\n                        let (cn, ck, cok)", "if false {\n                        let (cn, ck, cok)"),
    ("level-parse", F, "float_of(cn) / float_of(quantile.pow10(ck)), w, e);", "float_of(cn) / float_of(quantile.pow10(ck) + 1), w, e);"),
    ("corr-methods", F, "\"corr\") && cli.has(parsed, table, \"method\") && !bytes.equal(cli.text(args, parsed, table, \"method\"), \"pearson\") &&", "\"corr\") && false &&"),
    ("mean-guard", W, "    if n < 2 {\n        reason = \"too-few-rows\";\n    } else if mc != 0", "    if n < 1 {\n        reason = \"too-few-rows\";\n    } else if mc != 0"),
    ("mean-variance-zero", W, "} else if status == 2 {\n        ok = true;\n    } else {\n        let (q, qcode) = ttest.t_critical(", "} else if status == 1 {\n        ok = true;\n    } else {\n        let (q, qcode) = ttest.t_critical("),
    ("mean-df", W, "ttest.t_critical(lnum, lk, n - 1)", "ttest.t_critical(lnum, lk, n)"),
    ("mean-half", W, "half = q * float_of_bits(sb);", "half = q * float_of_bits(mb);"),
    ("mean-of", W, "put_mean_interval(heap, w, \"mean\", level", "put_mean_interval(heap, w, \"difference of means\", level"),
    ("mean-paired-of", W, "put_mean_interval(heap, w, \"difference of means\", level, float_of_bits(mb), half);", "put_mean_interval(heap, w, \"mean\", level, float_of_bits(mb), half);"),
    ("mean-estimate", W, "    if n >= 1 && mc == 0 {", "    if n >= 2 && mc == 0 {"),
    ("mean-pair-key", W, "w = json.put_key(heap, w, \"pair\");", "w = json.put_key(heap, w, \"pairs\");"),
    ("mean-missing", W, "w = json.put_int(heap, w, n);\n    w = json.put_key(heap, w, \"missing\");\n    w = json.put_int(heap, w, rows - n);", "w = json.put_int(heap, w, n);\n    w = json.put_key(heap, w, \"missing\");\n    w = json.put_int(heap, w, rows);"),
    ("diff-pooled-df", W, "df = float_of(r.n1 + r.n2 - 2);", "df = float_of(r.n1 + r.n2 - 1);"),
    ("diff-welch-df", W, "} else if r.status == 0 && r.fc == 0 {", "} else if r.status == 0 && r.fc == 1 {"),
    ("diff-variance-zero", W, "} else if r.status == 2 {\n        ok = true;", "} else if r.status == 1 {\n        ok = true;"),
    ("diff-half", W, "half = q * float_of_bits(r.sb);", "half = q * float_of_bits(r.eb);"),
    ("diff-estimate", W, "if r.status != 1 && r.ec == 0 {\n        w = json.put_float(heap, w, float_of_bits(r.eb));", "if r.status != 2 && r.ec == 0 {\n        w = json.put_float(heap, w, float_of_bits(r.eb));"),
    ("diff-n-a", W, "w = json.put_key(heap, w, \"a\");\n    w = json.put_int(heap, w, r.n1);", "w = json.put_key(heap, w, \"a\");\n    w = json.put_int(heap, w, r.n2);"),
    ("main-mean-dispatch", F, "if is_ci && bytes.equal(cli.operand(args, parsed, 1), \"mean\") {", "if is_ci && bytes.equal(cli.operand(args, parsed, 1), \"meanx\") {"),
    ("main-diff-dispatch", F, "} else if is_ci && bytes.equal(cli.operand(args, parsed, 1), \"diff\") {", "} else if is_ci && bytes.equal(cli.operand(args, parsed, 1), \"diffx\") {"),
    ("main-mean-pairs", F, "cli.has(parsed, table, \"pairs\"), false, true, w, e);", "cli.has(parsed, table, \"pairs\"), true, true, w, e);"),
    ("main-result-key", F, "rkey = \"intervals\";\n    }\n    var tname", "rkey = \"tests\";\n    }\n    var tname"),
    ("main-header-key", F, "hkey = \"ci\";\n        rkey = \"intervals\";\n    }\n    var tname", "hkey = \"test\";\n        rkey = \"intervals\";\n    }\n    var tname"),
    ("main-mean-name", F, "if ci {\n                                                    w = json.put_string(heap, w, \"mean\");", "if ci {\n                                                    w = json.put_string(heap, w, \"t1\");"),
    ("main-diff-name", F, "tname = \"diff\";", "tname = \"t2\";"),
    ("main-mu0", F, "if !paired && !ci {", "if !paired {"),
    ("main-mean-flags", F, "bytes.equal(cli.operand(args, parsed, 1), \"mean\") && !cli.has(parsed, table, \"cols\") && !cli.has(parsed, table, \"pairs\") {", "bytes.equal(cli.operand(args, parsed, 1), \"mean\") && !cli.has(parsed, table, \"cols\") || !cli.has(parsed, table, \"pairs\") {"),
    ("main-diff-flags", F, "bytes.equal(cli.operand(args, parsed, 1), \"diff\") && (!cli.has(parsed, table, \"cols\") || !cli.has(parsed, table, \"group\")) {", "bytes.equal(cli.operand(args, parsed, 1), \"diff\") && (!cli.has(parsed, table, \"cols\") && !cli.has(parsed, table, \"group\")) {"),
    ("corr-atanh-negative", W, "return 0.0 - 0.5 * math.log1p(2.0 * (0.0 - r) / (1.0 + r));", "return 0.0 - 0.5 * math.log1p(2.0 * (0.0 - r) / (1.0 - r));"),
    ("corr-atanh-positive", W, "return 0.5 * math.log1p(2.0 * r / (1.0 - r));", "return 0.5 * math.log1p(2.0 * r / (1.0 + r));"),
    ("corr-tanh-cut", W, "if a < 20.0 {", "if a < 10.0 {"),
    ("corr-tanh-sign", W, "if x < 0.0 {\n        return 0.0 - t;", "if x > 0.0 {\n        return 0.0 - t;"),
    ("corr-se", W, "let h = q / sqrt(float_of(n - 3));", "let h = q / sqrt(float_of(n - 2));"),
    ("corr-rows", W, "} else if n < 4 {", "} else if n < 3 {"),
    ("corr-perfect", W, "} else if exact == 1 {\n        reason = \"perfect-correlation\";", "} else if exact == 2 {\n        reason = \"perfect-correlation\";"),
    ("corr-constant", W, "if status != 0 {\n        reason = \"constant-column\";", "if status != 1 {\n        reason = \"constant-column\";"),
    ("corr-estimate", W, "    if status == 0 {\n        w = json.put_float(heap, w, r);", "    if status != 0 {\n        w = json.put_float(heap, w, r);"),
    ("corr-low", W, "lo = tanh(z - h);", "lo = tanh(z + h);"),
    ("corr-high", W, "hi = tanh(z + h);", "hi = tanh(z - h);"),
    ("corr-method", W, "\"correlation\", \"fisher-z\"", "\"correlation\", \"t\""),
    ("main-corr-dispatch", F, "} else if is_ci && bytes.equal(cli.operand(args, parsed, 1), \"corr\") {", "} else if is_ci && bytes.equal(cli.operand(args, parsed, 1), \"corrx\") {"),
    ("main-corr-header", F, "if ci {\n                                                w = json.put_key(heap, w, \"ci\");", "if !ci {\n                                                w = json.put_key(heap, w, \"ci\");"),
    ("main-corr-key", F, "w = json.put_key(heap, w, \"intervals\");\n                                            } else {", "w = json.put_key(heap, w, \"pairs\");\n                                            } else {"),
    ("main-corr-adjustment", F, "if !ci {\n                                                w = json.put_key(heap, w, \"adjustment\");", "if ci {\n                                                w = json.put_key(heap, w, \"adjustment\");"),
    ("main-corr-method", F, "\"corr\") && cli.has(parsed, table, \"method\") && !bytes.equal(cli.text(args, parsed, table, \"method\"), \"pearson\") {\n        e = rule.refuse(heap, e, \"args.bad-value\", \"ci corr", "\"corr\") && cli.has(parsed, table, \"method\") && false {\n        e = rule.refuse(heap, e, \"args.bad-value\", \"ci corr"),
    ("main-corr-level", F, "w = json.put_key(heap, w, \"level\");\n                                                w = json.put_float(heap, w, level);\n                                            }\n                                            w = json.put_key(heap, w, \"method\");", "w = json.put_key(heap, w, \"level\");\n                                                w = json.put_float(heap, w, 0.95);\n                                            }\n                                            w = json.put_key(heap, w, \"method\");"),
]

EXPLAINED = {
    "mean-variance-zero": "equivalent: with every value equal the standard error is 0, so the half-width q * 0 is 0 whether the branch that sets it is taken or the one that computes it (the quantile is found: the level is valid)",
    "cp-code": "unreachable guard (an argument, not a measurement): the bounds of section 9.2 of docs/numerics.md hold on every input of the sweeps, and no input of the gates or of the fuzz has made `beta_quantile` or `beta_quantile_upper` answer a nonzero code for shapes of 1 or more; the guard stays so that a failure is a reason, never a wrong endpoint",
}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_ci_prop.py"), "--bin", binary, "--files", str(files)], 1500),
                      ([sys.executable, str(ROOT / "scripts" / "gate_ci_mean.py"), "--bin", binary, "--files", str(files)], 1500),
                      ([sys.executable, str(ROOT / "scripts" / "gate_ci_corr.py"), "--bin", binary, "--files", str(files)], 1500),
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
        print("mutants_ci: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 100):
        print("mutants_ci: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_ci: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
