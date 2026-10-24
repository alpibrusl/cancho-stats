#!/usr/bin/env python3
"""Mutants of `stats test` (docs/test.md T-J): lib/stest/ttest.cho, the root of lib/acc/central.cho that its standard error needs, and the test part of tools/stats/main.cho. Each is applied to a copy, the
program rebuilt, scripts/gate_test.py and the refusals gate rerun; killed = a gate fails or the program traps. A mutant that does not compile is INVALID. The gates first run on the UNMUTATED build and the
script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.   python3 scripts/mutants_test.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
T = "lib/stest/ttest.cho"
C = "lib/acc/central.cho"
F = "tools/stats/main.cho"
A = "lib/stest/answer.cho"  # the pieces of the answers (put_pair, put_pair_kendall, put_t1), moved out of main.cho
M = [
    ("missing-cell-taken", T, "                    if key != rank.min_int() {\n                            let xb = bits_of(flt.float_of_key(key));\n                            psum.add(a[0..s1], 1, xb, wk);", "                    if key != rank.min_int() - 1 {\n                            let xb = bits_of(flt.float_of_key(key));\n                            psum.add(a[0..s1], 1, xb, wk);"),
    ("count-twice", T, "psum.add(a[s1..s1 + s2], 2, xb, wk);\n                            n = n + 1;\n                        }\n                    } else {", "psum.add(a[s1..s1 + s2], 2, xb, wk);\n                            n = n + 2;\n                        }\n                    } else {"),
    ("squares-as-first-power", T, "psum.add(a[s1..s1 + s2], 2, xb, wk);\n                            n = n + 1;\n                        }\n                    } else {", "psum.add(a[s1..s1 + s2], 1, xb, wk);\n                            n = n + 1;\n                        }\n                    } else {"),
    ("mu0-sign", T, "central.put(p, q, t[r_t..r_t + c], l1 + 3, smu >= 0);", "central.put(p, q, t[r_t..r_t + c], l1 + 3, smu < 0);"),
    ("mu0-times-n", T, "bigint.mul(t[r_t..r_t + c], t[r_mu..r_mu + c], l1, t[r_f..r_f + 3], 3);", "bigint.mul(t[r_t..r_t + c], t[r_mu..r_mu + c], l1, t[r_f + 16..r_f + 19], 3);"),
    ("mean-sign", T, "0 - 1074, sg1 < 0);", "0 - 1074, sg1 > 0);"),
    ("d-s1-sign", T, "central.put(p, q, t[r_s1..r_s1 + c], l1, sg1 < 0);", "central.put(p, q, t[r_s1..r_s1 + c], l1, sg1 > 0);"),
    ("mean-unit", T, "0 - 1074, sg1 < 0);", "0 - 1073, sg1 < 0);"),
    ("se-unit", T, "t[r_v..r_v + c], 6, 0 - 1074);", "t[r_v..r_v + c], 6, 0 - 1075);"),
    ("se-n-squared", T, "central.limbs_of(n * n, t[r_f + 8..r_f + 16]);", "central.limbs_of(n, t[r_f + 8..r_f + 16]);"),
    ("se-n-minus-one", T, "bigint.mul(t[r_v..r_v + c], t[r_f + 8..r_f + 11], 3, t[r_f + 16..r_f + 19], 3);", "bigint.mul(t[r_v..r_v + c], t[r_f + 8..r_f + 11], 3, t[r_f + 8..r_f + 11], 3);"),
    ("l-square-sign", T, "central.put(p, q, t[r_t..r_t + c], 2 * l1, true);", "central.put(p, q, t[r_t..r_t + c], 2 * l1, false);"),
    ("t-factor-n-minus-one", T, "bigint.mul(t[r_nu..r_nu + c], t[r_t..r_t + c], 2 * ld, t[r_f + 16..r_f + 19], 3);", "bigint.mul(t[r_nu..r_nu + c], t[r_t..r_t + c], 2 * ld, t[r_f..r_f + 3], 3);"),
    ("d-denominator-n", T, "bigint.mul(t[r_de..r_de + c], t[r_l..r_l + c], ll, t[r_f..r_f + 3], 3);", "bigint.mul(t[r_de..r_de + c], t[r_l..r_l + c], ll, t[r_f + 16..r_f + 19], 3);"),
    ("t-sign", T, "if sd < 0 {\n                                        if c1 == 0 {", "if sd > 0 {\n                                        if c1 == 0 {"),
    ("one-value", T, "if n < 2 {\n                                status = 1;", "if n < 1 {\n                                status = 1;"),
    ("constant-never", T, "if bigint.is_zero(t[r_l..r_l + c], c) {", "if bigint.is_zero(t[r_l..r_l + c], 0) {"),
    ("alpha-not-halved", T, "float_of(scale - num) / float_of(2 * scale)", "float_of(scale - num) / float_of(scale)"),
    ("root-exponent-ignored", C, "bigint.round_bits(contents(sw), sl, e - k, sticky)", "bigint.round_bits(contents(sw), sl, 0 - k, sticky)"),
    ("interval-width", A, "half = q * float_of_bits(sb);", "half = q * float_of_bits(sb) * 2.0;"),
    ("interval-df", A, "ttest.t_critical(lnum, lk, n - 1)", "ttest.t_critical(lnum, lk, n)"),
    ("p-df", A, "p_two_sided(float_of_bits(tb & 9223372036854775807), n - 1)", "p_two_sided(float_of_bits(tb & 9223372036854775807), n)"),
    ("p-uses-signed-t", A, "p_two_sided(float_of_bits(tb & 9223372036854775807), n - 1)", "p_two_sided(float_of_bits(tb), n - 1)"),
    ("effect-name", A, "w = json.put_string(heap, w, \"cohen_d\");", "w = json.put_string(heap, w, \"cohen_g\");"),
    ("method-name", A, "w = json.put_string(heap, w, \"one-sample t\");", "w = json.put_string(heap, w, \"t test\");"),
    ("df-printed", A, "    if n >= 2 {\n        w = json.put_int(heap, w, n - 1);", "    if n >= 2 {\n        w = json.put_int(heap, w, n);"),
    ("variance-zero-reason", A, "} else if status == 2 {\n        reason = \"variance-zero\";\n    } else if tc == 0 {", "} else if status == 2 {\n        reason = \"too-few-rows\";\n    } else if tc == 0 {"),
    ("missing-count", A, "w = json.put_int(heap, w, rows - n);", "w = json.put_int(heap, w, rows);"),
    ("top-level-echo", F, "w = json.put_float(heap, w, mu0);", "w = json.put_float(heap, w, 0.0);"),
    ("level-echo", F, "w = json.put_float(heap, w, level);\n                                                w = json.put_key(heap, w, \"adjustment\");", "w = json.put_float(heap, w, 0.95);\n                                                w = json.put_key(heap, w, \"adjustment\");"),
    ("cols-required", F, "} else if !cli.has(parsed, table, \"cols\") {\n            e = rule.refuse(heap, e, \"stats.test-columns\"", "} else if false {\n            e = rule.refuse(heap, e, \"stats.test-columns\""),
    ("kind-name", F, "!bytes.equal(cli.operand(args, parsed, 1), \"t1\") && !bytes", "!bytes.equal(cli.operand(args, parsed, 1), \"tX\") && !bytes"),
    ("level-one-allowed", F, "lv_num >= quantile.pow10(lv_k)", "lv_num > quantile.pow10(lv_k)"),
    ("level-zero-allowed", F, "!lv_ok || lv_num == 0 ||", "!lv_ok || lv_num == 99 ||"),
    ("mu0-check", F, "if mu_status != 0 {", "if mu_status == 99 {"),
    ("non-numeric-allowed", F, "if sel[j] == 1 && !numeric && ok {", "if sel[j] == 7 && !numeric && ok {"),
    ("dispatch", F, "                    } else if is_test {\n                        let (w2, e2) = test_one(", "                    } else if false {\n                        let (w2, e2) = test_one("),
    ("file-operand-index", F, "cli.operand_index(parsed, fidx)", "cli.operand_index(parsed, 1)"),
]

EXPLAINED = {
    "one-value": "equivalent: with one value L = 0, so the engine answers status 2 with every bit zero instead of status 1, and the answer is decided by n in put_t1 (statistic, df, p-value, effect size and interval are null, the reason is too-few-rows) before the status is read",
}


def control():
    """both gates on the UNMUTATED build must pass, or every kill below would be a gate that cannot run"""
    b = str(ROOT / "build" / "stats")
    g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_test.py"), "--bin", b, "--files", "60"], capture_output=True, text=True, timeout=900)
    g2 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe_refusals.py"), "--bin", b, "--fuzz", "150"], capture_output=True, text=True, timeout=900)
    return g.returncode == 0 and g2.returncode == 0


def main():
    for n, f, o, w in M:
        t = (ROOT / f).read_text()
        if t.count(o) != 1 or o == w:
            print("mutant %s: pattern occurs %d times" % (n, t.count(o)))
            return 1
    only = None
    if "--only" in sys.argv:
        only = set(sys.argv[sys.argv.index("--only") + 1].split(","))
    if "--check" in sys.argv:
        print("mutants_test: %d mutants, every pattern applies once" % len(M)); return 0
    if not control():
        print("mutants_test: the unmutated build fails the gates: refusing to count kills")
        return 1
    print("control (unmutated) passes both gates")
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
                bad.append(n); print("%-28s INVALID: %s" % (n, r.stderr[-160:].replace("\n", " "))); continue
            try:
                g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_test.py"), "--bin", str(d / "build" / "stats"), "--files", "60"], capture_output=True, text=True, timeout=300)
                dead = g.returncode != 0
                if not dead:
                    g2 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe_refusals.py"), "--bin", str(d / "build" / "stats"), "--fuzz", "150"], capture_output=True, text=True, timeout=600)
                    dead = g2.returncode != 0
            except subprocess.TimeoutExpired:
                dead = True
        if dead: killed += 1; print("%-28s killed" % n)
        elif n in EXPLAINED: print("%-28s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-28s SURVIVED" % n)
    total = len(M) if only is None else len([1 for n, f, o, w in M if n in only])
    print("mutants_test: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
