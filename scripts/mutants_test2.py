#!/usr/bin/env python3
"""Mutants of `stats test t2` and `paired` (docs/test.md T-J): lib/stest/welch.cho, lib/stest/answer2.cho, the paired part of lib/stest/ttest.cho, lib/csvx/pick.cho (--pairs) and the t2 and paired parts of tools/stats/main.cho.
Each is applied to a copy, the program rebuilt, scripts/gate_test2.py (t2 and paired) and the refusals gate rerun; killed = a gate fails or the program traps or hangs. A mutant that does not compile is INVALID. The gates
first run on the UNMUTATED build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.   python3 scripts/mutants_test2.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
W = "lib/stest/welch.cho"
B = "lib/stest/answer2.cho"
T = "lib/stest/ttest.cho"
P = "lib/csvx/pick.cho"
F = "tools/stats/main.cho"
A = "lib/stest/answer.cho"
M = [
    # the sums of the groups
    ("fill-missing-taken", W, "if key != rank.min_int() {\n            let xb", "if key != rank.min_int() - 1 {\n            let xb"),
    ("fill-count-twice", W, "psum.add(acc[s1..s1 + s2], 2, xb, wk);\n            n = n + 1;", "psum.add(acc[s1..s1 + s2], 2, xb, wk);\n            n = n + 2;"),
    ("fill-squares-as-first-power", W, "psum.add(acc[s1..s1 + s2], 2, xb, wk);", "psum.add(acc[s1..s1 + s2], 1, xb, wk);"),
    ("group-b-from-group-a", W, "n2 = fill(keys_b, ncols, col, rows_b, a[g..2 * g], wk);", "n2 = fill(keys_a, ncols, col, rows_a, a[g..2 * g], wk);"),
    # the difference of the means P = n2 S1a - n1 S1b
    ("p-sign-a", W, "l1 + 3, sga < 0);", "l1 + 3, sga > 0);"),
    ("p-sign-b", W, "l1 + 3, sgb >= 0);", "l1 + 3, sgb < 0);"),
    ("p-weight-a", W, "mulf(t, 8, 0, l1, f, 1, 3);", "mulf(t, 8, 0, l1, f, 0, 3);"),
    ("p-weight-b", W, "mulf(t, 8, 1, l1, f, 0, 3);", "mulf(t, 8, 1, l1, f, 1, 3);"),
    ("estimate-unit", W, "0 - 1074, sp < 0);", "0 - 1073, sp < 0);"),
    ("estimate-sign", W, "0 - 1074, sp < 0);", "0 - 1074, sp > 0);"),
    ("estimate-denominator", W, "f[8 * slot()..8 * slot() + 6], 6, 0 - 1074, sp < 0)", "f[10 * slot()..10 * slot() + 3], 3, 0 - 1074, sp < 0)"),
    # A and B (n S2 - S1^2 of each group)
    ("l-square-sign", W, "central.put(p, q, t[rt * c..rt * c + c], 2 * l1, true);", "central.put(p, q, t[rt * c..rt * c + c], 2 * l1, false);"),
    ("l-group-b-uses-n1", W, "lval(t, f, p, q, 1, 3, 5, 8, 1);", "lval(t, f, p, q, 1, 3, 5, 8, 0);"),
    ("l-group-b-uses-s2a", W, "lval(t, f, p, q, 1, 3, 5, 8, 1);", "lval(t, f, p, q, 1, 2, 5, 8, 1);"),
    # Welch: Q, the standard error, t, the degrees of freedom
    ("q-factors-swapped", W, "mulf(t, 8, 4, la, f, 7, 9);\n                        mulf(t, 9, 5, lb, f, 6, 9);", "mulf(t, 8, 4, la, f, 6, 9);\n                        mulf(t, 9, 5, lb, f, 7, 9);"),
    ("se-welch-unit", W, "f[12 * slot()..12 * slot() + 18], 18, 0 - 1074);", "f[12 * slot()..12 * slot() + 18], 18, 0 - 1075);"),
    ("se-welch-denominator", W, "f[12 * slot()..12 * slot() + 18], 18, 0 - 1074);", "f[13 * slot()..13 * slot() + 18], 18, 0 - 1074);"),
    ("t-welch-factor", W, "mulf(t, 11, 10, lp2, f, 9, 6);", "mulf(t, 11, 10, lp2, f, 11, 3);"),
    ("t-welch-denominator", W, "t[11 * c..12 * c], lp2 + 6, t[7 * c..8 * c], lq);", "t[11 * c..12 * c], lp2 + 6, t[12 * c..13 * c], len_of(t, 12));"),
    ("df-alpha-beta", W, "mulf(t, 9, 8, 2 * la, f, 13, 21);", "mulf(t, 9, 8, 2 * la, f, 14, 21);"),
    ("df-second-term", W, "mulf(t, 10, 8, 2 * lb, f, 14, 21);", "mulf(t, 10, 8, 2 * lb, f, 13, 21);"),
    ("df-numerator-factor", W, "mulf(t, 9, 8, 2 * lq, f, 9, 6);", "mulf(t, 9, 8, 2 * lq, f, 11, 3);"),
    ("df-denominator-sum", W, "sumr(t, 11, 9, 10);", "sumr(t, 11, 9, 9);"),
    # pooled: W, nu, t, the standard error
    ("pooled-w-weights", W, "mulf(t, 8, 4, la, f, 1, 3);\n                            mulf(t, 9, 5, lb, f, 0, 3);", "mulf(t, 8, 4, la, f, 0, 3);\n                            mulf(t, 9, 5, lb, f, 1, 3);"),
    ("pooled-se-denominator", W, "f[17 * slot()..17 * slot() + 15], 15, 0 - 1074);", "f[16 * slot()..16 * slot() + 12], 12, 0 - 1074);"),
    ("pooled-se-unit", W, "f[17 * slot()..17 * slot() + 15], 15, 0 - 1074);", "f[17 * slot()..17 * slot() + 15], 15, 0 - 1073);"),
    ("pooled-t-factor", W, "mulf(t, 8, 12, lw, f, 10, 3);", "mulf(t, 8, 12, lw, f, 11, 3);"),
    ("pooled-nu", W, "central.limbs_of(n1 + n2 - 2, f[11 * slot()..11 * slot() + 3]);", "central.limbs_of(n1 + n2 - 1, f[11 * slot()..11 * slot() + 3]);"),
    ("pooled-sum-of-sizes", W, "central.limbs_of(n1 + n2, f[10 * slot()..10 * slot() + 3]);", "central.limbs_of(n1 + n2 - 1, f[10 * slot()..10 * slot() + 3]);"),
    # Cohen's d and the ratio of the variances
    ("d-denominator", W, "mulf(t, 13, 12, lw, f, 8, 6);", "mulf(t, 13, 12, lw, f, 10, 3);"),
    ("d-numerator", W, "mulf(t, 11, 10, lp2, f, 11, 3);", "mulf(t, 11, 10, lp2, f, 10, 3);"),
    ("ratio-weights", W, "mulf(t, 8, 4, la, f, 18, 6);", "mulf(t, 8, 4, la, f, 19, 6);"),
    ("ratio-inverted", W, "central.ratio_bits(heap, t[8 * c..9 * c], len_of(t, 8), t[9 * c..10 * c], len_of(t, 9), 0, false)", "central.ratio_bits(heap, t[9 * c..10 * c], len_of(t, 9), t[8 * c..9 * c], len_of(t, 8), 0, false)"),
    ("ratio-zero-denominator", W, "if !bigint.is_zero(t[5 * c..6 * c], c) {", "if !bigint.is_zero(t[4 * c..5 * c], c) {"),
    # statuses and signs
    ("too-few-either", W, "if n1 < 2 || n2 < 2 {", "if n1 < 2 && n2 < 2 {"),
    ("variance-zero-never", W, "if bigint.is_zero(t[7 * c..8 * c], c) {", "if bigint.is_zero(t[7 * c..8 * c], 0) {"),
    ("t-sign", W, "if sp < 0 {\n                                if tc == 0 {", "if sp > 0 {\n                                if tc == 0 {"),
    ("d-sign", W, "if dc == 0 {\n                                    db = db | 0 - 9223372036854775807 - 1;", "if dc == 1 {\n                                    db = db | 0 - 9223372036854775807 - 1;"),
    # the answer
    ("pooled-df-value", B, "df = float_of(r.n1 + r.n2 - 2);", "df = float_of(r.n1 + r.n2 - 1);"),
    ("pooled-df-printed", B, "w = json.put_int(heap, w, r.n1 + r.n2 - 2);", "w = json.put_int(heap, w, r.n1 + r.n2 - 1);"),
    ("welch-df-scaled", B, "df = float_of_bits(r.fb);", "df = float_of_bits(r.fb) * 0.999;"),
    ("interval-df", B, "ttest.t_critical_real(lnum, lk, df)", "ttest.t_critical_real(lnum, lk, df + 1.0)"),
    ("interval-width", B, "half = q * float_of_bits(r.sb);", "half = q * float_of_bits(r.sb) * 2.0;"),
    ("p-df", B, "r.tb & 9223372036854775807), df);", "r.tb & 9223372036854775807), df + 1.0);"),
    ("p-uses-signed-t", B, "answer.p_two_sided_real(float_of_bits(r.tb & 9223372036854775807), df)", "answer.p_two_sided_real(float_of_bits(r.tb), df)"),
    ("interval-of", B, "w = json.put_string(heap, w, \"difference of means\");", "w = json.put_string(heap, w, \"mean\");"),
    ("methods-swapped", B, "w = json.put_string(heap, w, \"pooled\");\n    } else {\n        w = json.put_string(heap, w, \"welch\");", "w = json.put_string(heap, w, \"welch\");\n    } else {\n        w = json.put_string(heap, w, \"pooled\");"),
    ("ratio-null-shown", B, "if r.vstate == 0 {", "if r.vstate == 1 {"),
    ("sizes-swapped", B, "w = json.put_key(heap, w, \"a\");\n    w = json.put_int(heap, w, r.n1);", "w = json.put_key(heap, w, \"a\");\n    w = json.put_int(heap, w, r.n2);"),
    ("variance-zero-reason", B, "} else if r.status == 2 {\n        reason = \"variance-zero\";", "} else if r.status == 2 {\n        reason = \"too-few-rows\";"),
    ("effect-name", B, "w = json.put_string(heap, w, \"cohen_d\");", "w = json.put_string(heap, w, \"cohen_g\");"),
    ("zero-width-interval", B, "if r.status == 2 {\n            interval_ok = true;", "if r.status == 3 {\n            interval_ok = true;"),
    # the paired sums (ttest.cho)
    ("paired-b-sign", T, "let nyb = yb ^ (0 - 9223372036854775807 - 1);", "let nyb = yb;"),
    ("paired-cross-once", T, "xsum.add(a[s1..s1 + s2], xb, nyb, xk);\n                            xsum.add(a[s1..s1 + s2], xb, nyb, xk);", "xsum.add(a[s1..s1 + s2], xb, nyb, xk);"),
    ("paired-cross-sign", T, "xsum.add(a[s1..s1 + s2], xb, nyb, xk);\n                            xsum.add(a[s1..s1 + s2], xb, nyb, xk);", "xsum.add(a[s1..s1 + s2], xb, yb, xk);\n                            xsum.add(a[s1..s1 + s2], xb, yb, xk);"),
    ("paired-b-square-missing", T, "psum.add(a[s1..s1 + s2], 2, yb, wk);", "psum.add(a[s1..s1 + s2], 2, xb, wk);"),
    ("paired-a-in-s1-missing", T, "psum.add(a[0..s1], 1, xb, wk);\n                            psum.add(a[0..s1], 1, nyb, wk);", "psum.add(a[0..s1], 1, nyb, wk);"),
    ("paired-either-missing", T, "if key != rank.min_int() && key2 != rank.min_int() {", "if key != rank.min_int() || key2 != rank.min_int() {"),
    ("paired-interval-name", A, "w = json.put_string(heap, w, \"difference of means\");", "w = json.put_string(heap, w, \"mean\");"),
    ("paired-method-name", A, "w = json.put_string(heap, w, \"paired t\");", "w = json.put_string(heap, w, \"one-sample t\");"),
    ("paired-pair-order", A, "w = put_name(heap, w, names, ends, col);\n        w = put_name(heap, w, names, ends, col2);", "w = put_name(heap, w, names, ends, col2);\n        w = put_name(heap, w, names, ends, col);"),
    # --pairs
    ("pairs-colon", P, "if int_of(tok[k]) == 58 {", "if int_of(tok[k]) == 59 {"),
    ("pairs-empty-first-name", P, "if colon < 1 || colon + 1 >= len(tok) {", "if colon < 0 || colon + 1 >= len(tok) {"),
    ("pairs-empty-second-name", P, "if colon < 1 || colon + 1 >= len(tok) {", "if colon < 1 || colon + 1 > len(tok) {"),
    ("pairs-swapped", P, "pa[i] = ia;\n        pb[i] = ib;", "pa[i] = ib;\n        pb[i] = ia;"),
    ("pairs-ambiguity", P, "if ma > 1 || mb > 1 {", "if ma > 1 {"),
    ("pairs-no-match", P, "if ma == 0 || mb == 0 {", "if ma == 0 {"),
    # the verb
    ("pairs-limit-edge", F, "npairs > cli.nat(args, parsed, table, \"max-pairs\") {", "npairs >= cli.nat(args, parsed, table, \"max-pairs\") {"),
    ("pairs-second-column-unchecked", F, "contents(sw0)[contents(pbw)[k]] = 1;", "contents(sw0)[contents(pbw)[k]] = 0;"),
    ("paired-output-order", F, "contents(par)[k], contents(pbr)[k], info.rows, 0,", "contents(pbr)[k], contents(par)[k], info.rows, 0,"),
    ("paired-name", F, "w = json.put_string(heap, w, \"paired\");", "w = json.put_string(heap, w, \"t1\");"),
    ("pairs-required", F, "if !cli.has(parsed, table, \"pairs\") {", "if false {"),
    ("group-required", F, "&& !cli.has(parsed, table, \"group\") {", "&& false {"),
    ("group-levels-two", F, "if clean && nlevels != 2 {", "if clean && nlevels != 3 {"),
    ("group-a-is-second", F, "let key_a = map.key_at(gr.m, contents(gr.order)[0]);", "let key_a = map.key_at(gr.m, contents(gr.order)[1]);"),
    ("group-b-is-first", F, "let key_b = map.key_at(gr.m, contents(gr.order)[1]);", "let key_b = map.key_at(gr.m, contents(gr.order)[0]);"),
    ("rows-without-nulls", F, "info_a.rows + info_b.rows + gr.nulls", "info_a.rows + info_b.rows"),
    ("variance-flag", F, "bytes.equal(cli.text(args, parsed, table, \"variance\"), \"pooled\")", "bytes.equal(cli.text(args, parsed, table, \"variance\"), \"welch\")"),
    ("numeric-in-either-group", F, "let numeric = num_a[j] + num_b[j] > 0 && bad_a[j] == 0 && bad_b[j] == 0;", "let numeric = num_a[j] + num_b[j] > 0 && bad_a[j] == 0;"),
    ("numeric-needs-a-number", F, "let numeric = num_a[j] + num_b[j] > 0 && bad_a[j] == 0 && bad_b[j] == 0;", "let numeric = bad_a[j] == 0 && bad_b[j] == 0;"),
    ("t2-dispatch", F, "if is_test && bytes.equal(cli.operand(args, parsed, 1), \"t2\") {", "if false {"),
    ("t2-group-test-name", F, "w = json.put_string(heap, w, \"t2\");", "w = json.put_string(heap, w, \"t1\");"),
    ("kinds-accepted", F, "&& !bytes.equal(cli.operand(args, parsed, 1), \"paired\") && !bytes", "&& false && !bytes"),
]

EXPLAINED = {}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_test2.py"), "--bin", binary, "--kind", "t2", "--files", str(files)], 400),
                      ([sys.executable, str(ROOT / "scripts" / "gate_test2.py"), "--bin", binary, "--kind", "paired", "--files", str(files)], 400),
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
        print("mutants_test2: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 60):
        print("mutants_test2: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_test2: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
