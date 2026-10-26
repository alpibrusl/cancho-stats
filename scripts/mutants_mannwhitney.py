#!/usr/bin/env python3
"""Mutants of `stats test mannwhitney` (docs/rank.md R-H): lib/stest/mannwhitney.cho (U, the tie term, z, the exact p), lib/stest/exact.cho (the counts), lib/stest/answer7.cho (the entries), the mannwhitney parts of
lib/csvx/run.cho and of tools/stats/main.cho. Each is applied to a copy, the program rebuilt, scripts/gate_mannwhitney.py and the refusals gate rerun; killed = a gate fails or the program traps or hangs. A mutant that does
not compile is INVALID. The gates first run on the UNMUTATED build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.
    python3 scripts/mutants_mannwhitney.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
W = "lib/stest/mannwhitney.cho"
X = "lib/stest/exact.cho"
E = "lib/stest/answer7.cho"
R = "lib/csvx/run.cho"
F = "tools/stats/main.cho"
M = [
    # U and the effect size
    ("u-twice", W, "twice_u = contents(dw)[0] - na * (na + 1);", "twice_u = contents(dw)[0] - na * na;"),
    ("u-group", W, "twice_u = contents(dw)[0] - na * (na + 1);", "twice_u = contents(dw)[1] - na * (na + 1);"),
    ("u-centre", W, "let mm = twice_u - na * nb;", "let mm = twice_u - na * na;"),
    ("rb-denominator", W, "rb = float_of(mm) / float_of(na * nb);", "rb = float_of(mm) / float_of(na + nb);"),
    ("sizes", W, "nb = contents(nw)[1];", "nb = contents(nw)[0];"),
    ("groups-needed", W, "if na == 0 || nb == 0 {", "if na == 0 && nb == 0 {"),
    # the tie term and z
    ("e-cube", W, "welch.smul(f, 2, 1, 6, 0, 3);\n                        central.limbs_of(sum_t", "welch.smul(f, 2, 1, 6, 1, 3);\n                        central.limbs_of(sum_t"),
    ("e-plus-t", W, "bigint.add_in(t[8 * c..9 * c], c, f[3 * welch.slot()..3 * welch.slot() + 3], 3);", "bigint.add_in(t[8 * c..9 * c], c, f[0..3], 3);"),
    ("e-minus-n", W, "bigint.add_in(t[9 * c..10 * c], c, f[0..3], 3);", "bigint.add_in(t[9 * c..10 * c], c, f[3 * welch.slot()..3 * welch.slot() + 3], 3);"),
    ("e-zero", W, "let zero_e = bigint.is_zero(t[8 * c..9 * c], c);", "let zero_e = bigint.is_zero(t[9 * c..10 * c], c);"),
    ("continuity", W, "if cont == 1 && am > 0 {", "if cont == 0 && am > 0 {"),
    ("continuity-zero", W, "if cont == 1 && am > 0 {", "if cont == 1 && am >= 0 {"),
    ("z-square", W, "welch.smul(f, 5, 4, 3, 4, 3);", "welch.smul(f, 5, 4, 3, 5, 3);"),
    ("z-factor", W, "central.limbs_of(3 * m * (m - 1), f[6 * welch.slot()..6 * welch.slot() + 3]);", "central.limbs_of(3 * m * m, f[6 * welch.slot()..6 * welch.slot() + 3]);"),
    ("z-sizes", W, "central.limbs_of(na * nb, f[8 * welch.slot()..8 * welch.slot() + 3]);", "central.limbs_of(na + nb, f[8 * welch.slot()..8 * welch.slot() + 3]);"),
    ("z-sign", W, "if mm < 0 && zc == 0 && zb != 0 {", "if mm > 0 && zc == 0 && zb != 0 {"),
    # the exact p
    ("mode-auto", W, "use_exact = possible(keys, m);", "use_exact = false;"),
    ("mode-yes", W, "} else if mode == 1 {\n                            use_exact = true;", "} else if mode == 2 {\n                            use_exact = true;"),
    ("possible-limit", W, "if m > 60 {", "if m > 61 {"),
    ("possible-ties", W, "if keys[i] == keys[j] {", "if keys[i] == keys[j] + 1 {"),
    ("p-rank-sum", W, "exact.rank_sum_counts(heap, na, m, contents(dw)[0] / 2);", "exact.rank_sum_counts(heap, na, m, contents(dw)[0] / 2 + 1);"),
    ("p-tail", W, "if ge < lo {", "if ge > lo {"),
    ("p-double", W, "central.limbs_of(2 * lo, f[9 * welch.slot()..9 * welch.slot() + 3]);", "central.limbs_of(lo, f[9 * welch.slot()..9 * welch.slot() + 3]);"),
    ("p-cap", W, "if 2 * lo >= total {", "if 2 * lo > total {"),
    ("p-total", W, "central.limbs_of(total, f[10 * welch.slot()..10 * welch.slot() + 3]);", "central.limbs_of(total + 1, f[10 * welch.slot()..10 * welch.slot() + 3]);"),
    ("counts-ranks", X, "while s >= r {\n                    t[k * width + s]", "while s > r {\n                    t[k * width + s]"),
    ("counts-groups", X, "var k = na;\n            while k >= 1 {", "var k = na - 1;\n            while k >= 1 {"),
    ("counts-update", X, "t[k * width + s] = t[k * width + s] + t[(k - 1) * width + s - r];", "t[k * width + s] = t[k * width + s] + t[k * width + s - r];"),
    ("counts-le", X, "if s <= r_obs {", "if s < r_obs {"),
    ("counts-ge", X, "if s >= r_obs {", "if s > r_obs {"),
    ("counts-total", X, "total = total + c;", "total = total + c + c;"),
    ("counts-start", X, "t[0] = 1;\n        var r = 1;\n        while r <= nn {", "t[0] = 2;\n        var r = 1;\n        while r <= nn {"),
    # the entries
    ("entry-status", E, "    if r.status == 1 {\n        reason = \"too-few-rows\";", "    if r.status == 2 {\n        reason = \"too-few-rows\";"),
    ("entry-exact", E, "if r.status != 1 && r.exact == 1 {", "if r.status != 1 && r.exact == 0 {"),
    ("entry-normal-p", E, "p = 2.0 * sf;", "p = sf;"),
    ("entry-normal-cap", E, "if p > 1.0 {\n                p = 1.0;\n            }\n        } else {\n            reason = \"p-not-computed\";\n        }\n    }\n    w = json.put_key(heap, w, \"statistic\");", "if p > 2.0 {\n                p = 1.0;\n            }\n        } else {\n            reason = \"p-not-computed\";\n        }\n    }\n    w = json.put_key(heap, w, \"statistic\");"),
    ("entry-z-abs", E, "if az < 0.0 {\n            az = 0.0 - az;\n        }", "if az < 0.0 {\n            az = az;\n        }"),
    ("entry-u", E, "w = json.put_float(heap, w, float_of(r.twice_u) / 2.0);\n    } else {\n        w = json.put_null(heap, w);\n    }\n    w = json.put_key(heap, w, \"u\");", "w = json.put_float(heap, w, float_of(r.twice_u));\n    } else {\n        w = json.put_null(heap, w);\n    }\n    w = json.put_key(heap, w, \"u\");"),
    ("entry-ub", E, "float_of(2 * r.na * r.nb - r.twice_u) / 2.0", "float_of(r.na * r.nb - r.twice_u) / 2.0"),
    ("entry-continuity", E, "w = json.put_bool(heap, w, continuity);", "w = json.put_bool(heap, w, !continuity);"),
    ("entry-effect", E, "w = json.put_string(heap, w, \"rank_biserial\");", "w = json.put_string(heap, w, \"r\");"),
    ("entry-n", E, "w = json.put_int(heap, w, r.na);", "w = json.put_int(heap, w, r.nb);"),
    ("entry-ties", E, "w = json.put_int(heap, w, r.ties);", "w = json.put_int(heap, w, r.na);"),
    ("entry-beyond", E, "} else if r.zc != 0 {", "} else if r.zc != 1 {"),
    # the driver
    ("run-levels", R, "if clean && kind == 2 && k != 2 {", "if clean && kind == 2 && k != 3 {"),
    ("run-mode-yes", R, "if bytes.equal(exact_text, \"yes\") {", "if bytes.equal(exact_text, \"no\") {"),
    ("run-mode-no", R, "} else if bytes.equal(exact_text, \"no\") {\n        mode = 2;", "} else if bytes.equal(exact_text, \"yes\") {\n        mode = 2;"),
    ("run-continuity", R, "if bytes.equal(cont_text, \"none\") {\n        cont = 0;", "if bytes.equal(cont_text, \"auto\") {\n        cont = 0;"),
    ("run-precheck", R, "if clean && kind == 2 && bytes.equal(cli.text(args, parsed, table, \"exact\"), \"yes\") {", "if clean && kind == 2 && bytes.equal(cli.text(args, parsed, table, \"exact\"), \"yesx\") {"),
    ("run-precheck-count", R, "vec.size(gkr3) / nsel);", "vec.size(gkr3));"),
    ("run-kind-name", R, "w = json.put_string(heap, w, \"mannwhitney\");", "w = json.put_string(heap, w, \"kruskal\");"),
    ("run-stride", R, "let key = vec.get(gk, r * nsel + si);\n                        if key != rank.min_int() {\n                            contents(pkw)[m] = key;\n                            contents(pgw)[m] = vec.get(gi, r);", "let key = vec.get(gk, r * nsel);\n                        if key != rank.min_int() {\n                            contents(pkw)[m] = key;\n                            contents(pgw)[m] = vec.get(gi, r);"),
    ("main-group-required", F, "\"mannwhitney\") && !cli.has(parsed, table, \"group\") {", "\"mannwhitney\") && false {"),
    ("main-dispatch", F, "} else if bytes.equal(cli.operand(args, parsed, 1), \"mannwhitney\") {\n                            kind_of_groups = 2;", "} else if bytes.equal(cli.operand(args, parsed, 1), \"mannwhitneyx\") {\n                            kind_of_groups = 2;"),
]

EXPLAINED = {
    "counts-start": "equivalent: every count is doubled, and so is the total; the quotient 2 min(le, ge) / total is the same (the counts stay below 2^63: 2 C(60, 30) = 2.4e17)",
    "p-cap": "equivalent: when 2 min(le, ge) = total the quotient is exactly 1.0, the value the cap writes",
    "entry-normal-cap": "unreachable (argument, with z = 0 measured): 2 norm_sf(|z|) exceeds 1 only if norm_sf(|z|) exceeds 0.5, and norm_sf is decreasing from 0.5 at 0 (the gate meets z = 0, a p-value of exactly 1); the clamp stays as a guard of the contract p <= 1",
}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_mannwhitney.py"), "--bin", binary, "--files", str(files)], 1500),
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
        print("mutants_mannwhitney: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 100):
        print("mutants_mannwhitney: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_mannwhitney: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
