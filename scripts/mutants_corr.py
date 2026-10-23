#!/usr/bin/env python3
"""Mutants of `stats corr` (docs/corr.md C-I): lib/corr/pearson.cho, lib/corr/ranks.cho, the root of lib/acc/central.cho that the statistic needs, and the corr part of tools/stats/main.cho. Each is applied to a copy, the
driver rebuilt, scripts/gate_corr.py and the refusals gate rerun; killed = a gate fails or the program traps. A mutant that does not compile is INVALID. The gates first run on the UNMUTATED build and
the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.   python3 scripts/mutants_corr.py [--check]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
P = "lib/corr/pearson.cho"
C = "lib/acc/central.cho"
R = "lib/corr/ranks.cho"
F = "tools/stats/main.cho"
M = [
    ("mask-ignored", P, "if mask[i] != 0 && ka != rank.min_int()", "if ka != rank.min_int()"),
    ("missing-x-ignored", P, "if mask[i] != 0 && ka != rank.min_int() && kb != rank.min_int() {", "if mask[i] != 0 && kb != rank.min_int() {"),
    ("missing-y-ignored", P, "if mask[i] != 0 && ka != rank.min_int() && kb != rank.min_int() {", "if mask[i] != 0 && ka != rank.min_int() {"),
    ("key-of-x-column", P, "let ka = vec.get(keys, i * ncols + ca);", "let ka = vec.get(keys, i * ncols + cb);"),
    ("row-stride", P, "let kb = vec.get(keys, i * ncols + cb);", "let kb = vec.get(keys, i * ncols + cb + 1);"),
    ("sy-gets-x", P, "psum.add(a[oy..oy + s1], 1, yb, wk);", "psum.add(a[oy..oy + s1], 1, xb, wk);"),
    ("sxx-gets-y", P, "psum.add(a[oxx..oxx + s2], 2, xb, wk);", "psum.add(a[oxx..oxx + s2], 2, yb, wk);"),
    ("syy-gets-x", P, "psum.add(a[oyy..oyy + s2], 2, yb, wk);", "psum.add(a[oyy..oyy + s2], 2, xb, wk);"),
    ("sxy-gets-yy", P, "xsum.add(a[oxy..oxy + s2], xb, yb, wk);", "xsum.add(a[oxy..oxy + s2], yb, yb, wk);"),
    ("count-twice", P, "                    n = n + 1;", "                    n = n + 2;"),
    ("lx-square-sign", P, "central.put(p, q, t[r_t..r_t + c], 2 * l1, true);\n                        central.settle(p, q, t[r_lx..r_lx + c]);", "central.put(p, q, t[r_t..r_t + c], 2 * l1, false);\n                        central.settle(p, q, t[r_lx..r_lx + c]);"),
    ("lx-uses-syy", P, "bigint.mul(t[r_t..r_t + c], t[r_sxx..r_sxx + c], l2, t[r_f..r_f + 8], 3);", "bigint.mul(t[r_t..r_t + c], t[r_syy..r_syy + c], l2, t[r_f..r_f + 8], 3);"),
    ("ly-square-sign", P, "central.put(p, q, t[r_t..r_t + c], 2 * l1, true);\n                        central.settle(p, q, t[r_ly..r_ly + c]);", "central.put(p, q, t[r_t..r_t + c], 2 * l1, false);\n                        central.settle(p, q, t[r_ly..r_ly + c]);"),
    ("constant-needs-both", P, "if bigint.is_zero(t[r_lx..r_lx + c], c) || bigint.is_zero(t[r_ly..r_ly + c], c) {", "if bigint.is_zero(t[r_lx..r_lx + c], c) && bigint.is_zero(t[r_ly..r_ly + c], c) {"),
    ("n-sxy-sign", P, "central.put(p, q, t[r_t..r_t + c], l2 + 3, sxy < 0);", "central.put(p, q, t[r_t..r_t + c], l2 + 3, sxy > 0);"),
    ("n-sxsy-sign", P, "(sx < 0) == (sy < 0));", "(sx < 0) != (sy < 0));"),
    ("exact-never", P, "if bigint.compare(t[r_nn..r_nn + c], lnn, t[r_d..r_d + c], ld) == 0 {", "if bigint.compare(t[r_nn..r_nn + c], lnn, t[r_d..r_d + c], ld) > 0 {"),
    ("exact-constant", P, "rb = 4607182418800017408;", "rb = 4607182418800017409;"),
    ("exact-sign", P, "if sn < 0 {\n                                    rb = rb | 0 - 9223372036854775807 - 1;", "if sn > 0 {\n                                    rb = rb | 0 - 9223372036854775807 - 1;"),
    ("r-sign", P, "if sn < 0 && c1 == 0 {", "if sn > 0 && c1 == 0 {"),
    ("t-needs-four-rows", P, "                                if n >= 3 {\n                                    // t^2", "                                if n >= 4 {\n                                    // t^2"),
    ("t-df", P, "central.limbs_of(n - 2, t[r_f..r_f + 8]);", "central.limbs_of(n - 1, t[r_f..r_f + 8]);"),
    ("t-denominator-sum", P, "bigint.sub_in(contents(denw), nd, t[r_nn..r_nn + c], lnn);", "bigint.add_in(contents(denw), nd, t[r_nn..r_nn + c], lnn);"),
    ("nn-length", P, "let lnn = 2 * ln;", "let lnn = ln;"),
    ("d-length", P, "let ld = lx + ly;", "let ld = lx;"),
    ("root-limbs-fixed", C, "    var sl = top / 26 + 2;", "    var sl = 8;"),
    ("statistic-sign", F, "if rb < 0 {\n            tbits = tb |", "if rb > 0 {\n            tbits = tb |"),
    ("df", F, "w = json.put_int(heap, w, n - 2);", "w = json.put_int(heap, w, n - 1);"),
    ("p-one-sided", F, "var p = 2.0 * sf;", "var p = sf;"),
    ("p-clamp", F, "if p > 1.0 {\n        p = 1.0;", "if p > 2.0 {\n        p = 1.0;"),
    ("two-rows-p", F, "if n < 3 {\n            reason = \"too-few-rows\";", "if n < 2 {\n            reason = \"too-few-rows\";"),
    ("exact-p-value", F, "w = json.put_float(heap, w, 0.0);", "w = json.put_float(heap, w, 1.0);"),
    ("beyond-double-p", F, "reason = \"statistic-beyond-double\";\n        } else {", "reason = \"statistic-beyond-double\";\n            have_p = true;\n        } else {"),
    ("listwise-flag", F, "let listwise = bytes.equal(cli.text(args, parsed, table, \"missing\"), \"listwise\");", "let listwise = false;"),
    ("listwise-mask", F, "contents(mw)[i] = 0;", "contents(mw)[i] = 1;"),
    ("pair-with-itself", F, "var ib = ia + 1;", "var ib = ia;"),
    ("text-column-taken", F, "let numeric = contents(c0.numeric)[j] > 0 && contents(c0.bad)[j] == 0;", "let numeric = contents(c0.numeric)[j] > 0;"),
    ("named-text-column", F, "} else if contents(sw1)[j] == 1 && !numeric && cols_ok {", "} else if contents(sw1)[j] == 7 && !numeric && cols_ok {"),
    ("one-column-is-enough", F, "if cols_ok && m < 2 {", "if cols_ok && m < 1 {"),
    ("pairs-limit-edge", F, "if cols_ok && pairs > cli.nat(args, parsed, table, \"max-pairs\") {", "if cols_ok && pairs >= cli.nat(args, parsed, table, \"max-pairs\") {"),
    ("pairs-ceiling", F, "cli.nat(args, parsed, table, \"max-pairs\") > 32640", "cli.nat(args, parsed, table, \"max-pairs\") > 32641"),
    ("pairs-floor", F, "cli.nat(args, parsed, table, \"max-pairs\") < 1 {", "cli.nat(args, parsed, table, \"max-pairs\") < 0 {"),
    ("verb-corr", F, "&& !bytes.equal(cli.operand(args, parsed, 0), \"corr\") {", "&& !bytes.equal(cli.operand(args, parsed, 0), \"corr \") {"),
    ("rank-plus-one", R, "float_of(lo + hi + 1)", "float_of(lo + hi)"),
    ("rank-upper-is-lower", R, "let hi = upper(contents(sw), m, key);", "let hi = lower(contents(sw), m, key);"),
    ("upper-strict", R, "if a[mid] <= key {", "if a[mid] < key {"),
    ("lower-nonstrict", R, "if a[mid] < key {", "if a[mid] <= key {"),
    ("used-ignores-mask", R, "if mask[i] != 0 && vec.get(keys, i * ncols + ca) != rank.min_int() && vec.get(keys, i * ncols + cb) != rank.min_int() {\n            used = used + 1;", "if vec.get(keys, i * ncols + ca) != rank.min_int() && vec.get(keys, i * ncols + cb) != rank.min_int() {\n            used = used + 1;"),
    ("gather-ignores-partner", R, "if mask[i] != 0 && vec.get(keys, i * ncols + ca) != rank.min_int() && vec.get(keys, i * ncols + cb) != rank.min_int() {\n                                contents(gw)[m] =", "if mask[i] != 0 && vec.get(keys, i * ncols + ca) != rank.min_int() {\n                                contents(gw)[m] ="),
    ("tie-group-size", R, "if end - at > 1 {", "if end - at > 2 {"),
    ("tie-run-test", R, "contents(sw)[end] == contents(sw)[at]", "contents(sw)[end] != contents(sw)[at]"),
    ("rank-column-b", R, "col = cb;", "col = ca;"),
    ("ties-swapped", R, "if side == 0 {\n            ties0 = groups;", "if side == 1 {\n            ties0 = groups;"),
    ("table-side", R, "vec.set(tbl, 2 * i + side,", "vec.set(tbl, 2 * i,"),
    ("not-sorted", R, "rank.sort_range(contents(sw), contents(tw), contents(cw), 0, m);", "rank.sort_range(contents(sw), contents(tw), contents(cw), 0, 0);"),
    ("spearman-flag", F, "\"method\"), \"spearman\") {", "\"method\"), \"spearmen\") {"),
    ("rho-key", F, "w = json.put_key(heap, w, \"rho\");", "w = json.put_key(heap, w, \"r\");"),
    ("p-method-spearman", F, "w = json.put_string(heap, w, \"t-approximation\");", "w = json.put_string(heap, w, \"t\");"),
    ("all-tied-reason", F, "reason = \"all-tied\";", "reason = \"constant-column\";"),
    ("ties-y", F, "w = json.put_int(heap, w, ties_y);", "w = json.put_int(heap, w, ties_x);"),
    ("ties-rule", F, "w = json.put_string(heap, w, \"average\");", "w = json.put_string(heap, w, \"ordinal\");"),
    ("pearson-on-store", F, "pearson.pair(heap, tr, 2, 0, 1, rows, contents(or));", "pearson.pair(heap, store, ncols, ja, jb, rows, mask);"),
]

EXPLAINED = {
    "rank-plus-one": "equivalent: lo + hi + 1 is twice the average rank, and lo + hi is twice the average rank minus one for every row; Pearson's coefficient and its statistic do not change when a constant is added to all the values of a column",
    "used-ignores-mask": "equivalent: `used` only sizes the three buffers of the column, the gather loop below counts the rows it takes in `m`, so counting the rows the mask excludes allocates more than is used and changes no value",
    "p-clamp": "defensive: 2 * sf exceeds 1 only if the Student t survival function answers more than 1/2 for a positive t, which no gate case reached (the gate checks every p-value is in [0, 1] and has a case with r = 0 and one with r about 1e-17); a p-value above 1 by one rounding would also lie inside the 1e-12 bound of the p-value gate, so the clamp cannot be seen by it",
}


def control():
    """both gates on the UNMUTATED build must pass, or every kill below would be a gate that cannot run"""
    b = str(ROOT / "build" / "stats")
    g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_corr.py"), "--bin", b, "--files", "60"], capture_output=True, text=True, timeout=900)
    g1 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_corr.py"), "--bin", b, "--files", "60", "--method", "spearman"], capture_output=True, text=True, timeout=900)
    g2 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe_refusals.py"), "--bin", b, "--fuzz", "150"], capture_output=True, text=True, timeout=900)
    return g.returncode == 0 and g1.returncode == 0 and g2.returncode == 0


def main():
    for n, f, o, w in M:
        t = (ROOT / f).read_text()
        if t.count(o) != 1 or o == w:
            print("mutant %s: pattern occurs %d times" % (n, t.count(o)))
            return 1
    if "--check" in sys.argv:
        print("mutants_corr: %d mutants, every pattern applies once" % len(M)); return 0
    if not control():
        print("mutants_corr: the unmutated build fails the gates: refusing to count kills")
        return 1
    print("control (unmutated) passes both gates")
    bad, killed = [], 0
    for n, f, o, w in M:
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
                g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_corr.py"), "--bin", str(d / "build" / "stats"), "--files", "60"], capture_output=True, text=True, timeout=300)
                dead = g.returncode != 0
                if not dead:
                    g1 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_corr.py"), "--bin", str(d / "build" / "stats"), "--files", "60", "--method", "spearman"], capture_output=True, text=True, timeout=300)
                    dead = g1.returncode != 0
                if not dead:
                    g2 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe_refusals.py"), "--bin", str(d / "build" / "stats"), "--fuzz", "150"], capture_output=True, text=True, timeout=600)
                    dead = g2.returncode != 0
            except subprocess.TimeoutExpired:
                dead = True
        if dead: killed += 1; print("%-28s killed" % n)
        elif n in EXPLAINED: print("%-28s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-28s SURVIVED" % n)
    print("mutants_corr: %d killed of %d, %d survived and explained, %d unexplained" % (killed, len(M), len(M) - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
