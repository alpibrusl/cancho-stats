#!/usr/bin/env python3
"""Mutants of `stats test wilcoxon` (docs/rank.md R-H): lib/stest/wilcoxon.cho (the exact differences, the ranks, the three methods for zeros, the normal score, the exact p), lib/stest/exact.cho (the counts of signed-rank
sums), lib/stest/answer8.cho (the entries) and the wilcoxon parts of tools/stats/main.cho. Each is applied to a copy, the program rebuilt, scripts/gate_wilcoxon.py and the refusals gate rerun; killed = a gate fails or the
program traps or hangs. A mutant that does not compile is INVALID. The gates first run on the UNMUTATED build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.
    python3 scripts/mutants_wilcoxon.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
W = "lib/stest/wilcoxon.cho"
X = "lib/stest/exact.cho"
E = "lib/stest/answer8.cho"
F = "tools/stats/main.cho"
M = [
    # the exact differences and their order
    ("diff-sign", W, "dd.two_sum(flt.float_of_key(ka), 0.0 - b);", "dd.two_sum(flt.float_of_key(ka), b);"),
    ("diff-constant", W, "var b = mu0;", "var b = 0.0;"),
    ("diff-overflow", W, "if s0 > 1.7976931348623157e308 || s0 < 0.0 - 1.7976931348623157e308 {", "if s0 > 1.7976931348623157e308 {"),
    ("diff-present", W, "if ka != rank.min_int() && kb != rank.min_int() {", "if ka != rank.min_int() || kb != rank.min_int() {"),
    ("diff-error-key", W, "a2[m] = rank.key_of(e);", "a2[m] = 0;"),
    ("diff-error-sign", W, "a2[m] = rank.key_of(0.0 - e);", "a2[m] = rank.key_of(e);"),
    ("diff-magnitude-sign", W, "a1[m] = rank.key_of(0.0 - s);", "a1[m] = rank.key_of(s);"),
    ("diff-sign-flag", W, "sgn[m] = 0 - 1;", "sgn[m] = 1;"),
    ("order-second-key", W, "return k2[x] < k2[y];", "return k2[x] > k2[y];"),
    ("order-first-key", W, "return k1[x] < k1[y];", "return k1[x] > k1[y];"),
    ("order-tie-equal", W, "a1[ix[end]] == a1[ix[at]] && a2[ix[end]] == a2[ix[at]]", "a1[ix[end]] == a1[ix[at]]"),
    ("order-stable", W, "if less(k1, k2, idx[y], idx[x]) {", "if !less(k1, k2, idx[x], idx[y]) {"),
    # the zeros
    ("zeros-ranked", W, "if zmode >= 1 {\n                                                a1[m] = 0;", "if zmode >= 2 {\n                                                a1[m] = 0;"),
    ("zeros-zsplit-sums", W, "dp = dp + zeros * (zeros + 1) / 2;\n                                        dm = dm + zeros * (zeros + 1) / 2;", "dp = dp + zeros * (zeros + 1) / 2;\n                                        dm = dm;"),
    ("zeros-zsplit-tie", W, "if sgn[ix[at]] == 0 && zmode != 2 {", "if sgn[ix[at]] == 0 && zmode != 1 {"),
    ("zeros-normal-size", W, "if zmode >= 1 {\n                                        nn = used;", "if zmode >= 2 {\n                                        nn = used;"),
    ("zeros-pratt-mean", W, "mean2 = (nn * (nn + 1) - zeros * (zeros + 1)) / 2;", "mean2 = nn * (nn + 1) / 2;"),
    ("zeros-pratt-variance", W, "welch.smul(f, 12, 11, 6, 10, 3);", "welch.smul(f, 12, 11, 6, 9, 3);"),
    ("zeros-only", W, "if nz == 0 && zmode == 0 {", "if nz == 0 && zmode == 2 {"),
    # ranks, sums, tie term
    ("rank-doubled", W, "let rk = at + end + 1;", "let rk = at + end;"),
    ("rank-sum-plus", W, "if sgn[ix[q]] > 0 {\n                                                dp = dp + rk;", "if sgn[ix[q]] < 0 {\n                                                dp = dp + rk;"),
    ("tie-size", W, "var counts_tie = size > 1;", "var counts_tie = size > 2;"),
    ("tie-cube", W, "welch.smul(f, 2, 1, 6, 0, 3);\n                                            bigint.add_in(t[6 * c", "welch.smul(f, 2, 1, 6, 1, 3);\n                                            bigint.add_in(t[6 * c"),
    ("tie-sum", W, "sum_t = sum_t + size;", "sum_t = sum_t + 1;"),
    ("rb-denominator", W, "rb = float_of(dp - dm) / float_of(dp + dm);", "rb = float_of(dp - dm) / float_of(dp);"),
    # the normal score
    ("mean-wilcox", W, "var mean2 = nn * (nn + 1) / 2;", "var mean2 = nn * (nn + 1) / 4;"),
    ("variance-factor", W, "welch.smul(f, 7, 6, 6, 5, 3);", "welch.smul(f, 7, 6, 6, 4, 3);"),
    ("variance-doubled", W, "bigint.add_in(t[c..2 * c], c, t[0..c], c);", "bigint.add_in(t[c..2 * c], c, t[7 * c..8 * c], c);"),
    ("variance-ties-plus", W, "bigint.add_in(t[c..2 * c], c, f[13 * welch.slot()..13 * welch.slot() + 3], 3);", "bigint.add_in(t[c..2 * c], c, f[3 * welch.slot()..3 * welch.slot() + 3], 3);"),
    ("variance-ties-minus", W, "bigint.sub_in(t[c..2 * c], c, t[6 * c..7 * c], c);", "bigint.sub_in(t[c..2 * c], c, t[0..c], c);"),
    ("variance-zero", W, "if status == 0 && bigint.is_zero(t[c..2 * c], c) {", "if status == 0 && bigint.is_zero(t[0..c], c) {"),
    ("continuity", W, "if cont == 1 && am > 0 {", "if cont == 0 && am > 0 {"),
    ("continuity-zero", W, "if cont == 1 && am > 0 {", "if cont == 1 && am >= 0 {"),
    ("z-twelve", W, "central.limbs_of(12, f[2 * welch.slot()..2 * welch.slot() + 3]);", "central.limbs_of(6, f[2 * welch.slot()..2 * welch.slot() + 3]);"),
    ("z-square", W, "welch.smul(f, 1, 0, 3, 0, 3);\n                                        central.limbs_of(12", "welch.smul(f, 1, 0, 3, 3, 3);\n                                        central.limbs_of(12"),
    ("z-sign", W, "if dp < mean2 && zc == 0 && zb != 0 {", "if dp > mean2 && zc == 0 && zb != 0 {"),
    # the exact p
    ("exact-eligible-zeros", W, "let eligible = (zeros == 0 || zmode == 0) && ties == 0 && nz <= 50 && nz > 0 && (zmode == 0 || zeros == 0);", "let eligible = ties == 0 && nz <= 50 && nz > 0;"),
    ("exact-limit", W, "nz <= 50 && nz > 0", "nz <= 51 && nz > 0"),
    ("exact-ties", W, "let eligible = (zeros == 0 || zmode == 0) && ties == 0", "let eligible = (zeros == 0 || zmode == 0)"),
    ("exact-mode", W, "} else if mode == 1 {\n                                        use_exact = true;", "} else if mode == 2 {\n                                        use_exact = true;"),
    ("exact-statistic", W, "exact.signed_rank_counts(heap, nz, dp / 2);", "exact.signed_rank_counts(heap, nz, dp / 2 + 1);"),
    ("exact-tail", W, "if ge < lo {", "if ge > lo {"),
    ("exact-total", W, "let total = 1 << nz;", "let total = 1 << (nz - 1);"),
    ("exact-double", W, "central.limbs_of(2 * lo, f[0..3]);", "central.limbs_of(lo, f[0..3]);"),
    ("exact-cap", W, "if 2 * lo >= total {", "if 2 * lo > total {"),
    ("counts-ranks", X, "while s >= r {\n                t[s] = t[s] + t[s - r];", "while s > r {\n                t[s] = t[s] + t[s - r];"),
    ("counts-update", X, "t[s] = t[s] + t[s - r];", "t[s] = t[s] + t[s - 1];"),
    ("counts-le", X, "if s <= w_obs {\n                le = le + t[s];", "if s < w_obs {\n                le = le + t[s];"),
    ("counts-ge", X, "if s >= w_obs {\n                ge = ge + t[s];", "if s > w_obs {\n                ge = ge + t[s];"),
    ("counts-start", X, "t[0] = 1;\n        var r = 1;\n        while r <= n {", "t[0] = 1;\n        var r = 2;\n        while r <= n {"),
    # the entries
    ("entry-status", E, "    if r.status == 1 {\n        reason = \"too-few-rows\";", "    if r.status == 2 {\n        reason = \"too-few-rows\";"),
    ("entry-overflow", E, "} else if r.status == 3 {\n        reason = \"difference-beyond-double\";", "} else if r.status == 2 {\n        reason = \"difference-beyond-double\";"),
    ("entry-answered", E, "let answered = r.status == 0 || r.status == 2;", "let answered = r.status == 0;"),
    ("entry-exact-p", E, "if r.status == 0 || r.status == 2 {\n        if r.exact == 1 {", "if r.status == 0 {\n        if r.exact == 1 {"),
    ("entry-normal-p", E, "p = 2.0 * sf;", "p = sf;"),
    ("entry-statistic", E, "w = json.put_float(heap, w, float_of(r.dp) / 2.0);\n    } else {\n        w = json.put_null(heap, w);\n    }\n    w = json.put_key(heap, w, \"w\");", "w = json.put_float(heap, w, float_of(r.dm) / 2.0);\n    } else {\n        w = json.put_null(heap, w);\n    }\n    w = json.put_key(heap, w, \"w\");"),
    ("entry-minus", E, "w = json.put_float(heap, w, float_of(r.dm) / 2.0);", "w = json.put_float(heap, w, float_of(r.dp) / 2.0);"),
    ("entry-effect-null", E, "if answered && r.dp + r.dm > 0 {", "if answered {"),
    ("entry-n-used", E, "w = json.put_int(heap, w, r.n);", "w = json.put_int(heap, w, r.nz);"),
    ("entry-n-zeros", E, "w = json.put_int(heap, w, r.zeros);", "w = json.put_int(heap, w, r.nz);"),
    ("entry-zero-method", E, "} else if zmode == 1 {\n        w = json.put_string(heap, w, \"pratt\");", "} else if zmode == 2 {\n        w = json.put_string(heap, w, \"pratt\");"),
    ("entry-continuity", E, "w = json.put_bool(heap, w, cont == 1);", "w = json.put_bool(heap, w, cont == 0);"),
    ("entry-pair", E, "w = answer.put_name(heap, w, names, ends, col2);", "w = answer.put_name(heap, w, names, ends, col);"),
    ("entry-exact-ok", E, "return r.status != 0 || r.exact == 1;", "return r.exact == 1;"),
    # the driver
    ("main-header", F, "if wil {\n                                                    w = json.put_string(heap, w, \"wilcoxon\");", "if !wil {\n                                                    w = json.put_string(heap, w, \"wilcoxon\");"),
    ("main-zeros-pratt", F, "if bytes.equal(cli.text(args, parsed, table, \"zeros\"), \"pratt\") {\n                        zmode = 1;", "if bytes.equal(cli.text(args, parsed, table, \"zeros\"), \"zsplit\") {\n                        zmode = 1;"),
    ("main-zeros-zsplit", F, "} else if bytes.equal(cli.text(args, parsed, table, \"zeros\"), \"zsplit\") {\n                        zmode = 2;", "} else if bytes.equal(cli.text(args, parsed, table, \"zeros\"), \"pratt\") {\n                        zmode = 2;"),
    ("main-exact-yes", F, "if bytes.equal(cli.text(args, parsed, table, \"exact\"), \"yes\") {\n                        xmode = 1;", "if bytes.equal(cli.text(args, parsed, table, \"exact\"), \"no\") {\n                        xmode = 1;"),
    ("main-exact-no", F, "} else if bytes.equal(cli.text(args, parsed, table, \"exact\"), \"no\") {\n                        xmode = 2;", "} else if bytes.equal(cli.text(args, parsed, table, \"exact\"), \"yes\") {\n                        xmode = 2;"),
    ("main-continuity", F, "if bytes.equal(cli.text(args, parsed, table, \"continuity\"), \"none\") {\n                        cmode = 0;", "if bytes.equal(cli.text(args, parsed, table, \"continuity\"), \"auto\") {\n                        cmode = 0;"),
    ("main-precheck", F, "if cols_ok && wil && xmode == 1 {", "if cols_ok && wil && xmode == 7 {"),
    ("main-precheck-pairs", F, "!answer8.exact_ok(heap, sb0, nc, contents(par0)[k], contents(pbr0)[k], info.rows, 0.0, zmode, cmode)", "!answer8.exact_ok(heap, sb0, nc, contents(par0)[k], contents(par0)[k], info.rows, 0.0, zmode, cmode)"),
    ("main-pairs-or-cols", F, "if cli.has(parsed, table, \"pairs\") == cli.has(parsed, table, \"cols\") {\n                e = rule.refuse(heap, e, \"stats.test-columns\", \"wilcoxon needs", "if cli.has(parsed, table, \"pairs\") && cli.has(parsed, table, \"cols\") {\n                e = rule.refuse(heap, e, \"stats.test-columns\", \"wilcoxon needs"),
    ("main-mu0", F, "answer8.put_wilcoxon(heap, w, buffer.bytes(nb), eb, sb, nc, j, 0 - 1, info.rows, mu0, zmode, xmode, cmode);", "answer8.put_wilcoxon(heap, w, buffer.bytes(nb), eb, sb, nc, j, 0 - 1, info.rows, 0.0, zmode, xmode, cmode);"),
    ("main-level-hidden", F, "                                                if !wil {\n                                                    w = json.put_key(heap, w, \"level\");", "                                                if wil {\n                                                    w = json.put_key(heap, w, \"level\");"),
]

EXPLAINED = {
    "variance-zero": "equivalent: Vnum = 2 S + sum t - sum t^3 is zero only when S is (S = 0 is the one case, `pratt` with only zero differences; with a nonzero difference Vnum is positive, as the variance of the signed ranks is)",
    "order-stable": "equivalent: equal keys are one tie group, ranked alike whichever way they were ordered among themselves",
    "exact-cap": "equivalent: when 2 min(le, ge) = total the quotient is exactly 1.0, the value the cap writes",
    "entry-exact-p": "equivalent: a variance of zero (status 2) comes only with no nonzero difference, and then no exact p-value is computed (nz > 0 is a condition of it)",
}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_wilcoxon.py"), "--bin", binary, "--files", str(files)], 1500),
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
        print("mutants_wilcoxon: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 100):
        print("mutants_wilcoxon: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_wilcoxon: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
