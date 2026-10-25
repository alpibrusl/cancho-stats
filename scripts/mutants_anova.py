#!/usr/bin/env python3
"""Mutants of `stats test anova` (docs/rank.md R-H): lib/stest/anova.cho (the exact sums of squares, the effect sizes, Welch's weights, distances and double-double core), lib/stest/answer5.cho (the entries) and lib/csvx/run.cho (the
reading of the groups, the checks). Each is applied to a copy, the program rebuilt, scripts/gate_anova.py and the refusals gate rerun; killed = a gate fails or the program traps or hangs. A mutant that does not compile is
INVALID. The gates first run on the UNMUTATED build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.   python3 scripts/mutants_anova.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
A = "lib/stest/anova.cho"
F = "tools/stats/main.cho"
W = "lib/stest/welchf.cho"
E = "lib/stest/answer5.cho"
R = "lib/csvx/run.cho"
M = [
    # Welch's statistic: the double-double conversion, the pairs, the tail
    ("ratio-extra-bits", W, "var k = 120 + bd - bn;", "var k = 20 + bd - bn;"),
    ("ratio-limb-scale", W, "h = h * 67108864.0;\n                l = l * 67108864.0;", "h = h * 33554432.0;\n                l = l * 33554432.0;"),
    ("ratio-cut-bits", W, "let sh = bn - bd - 140;", "let sh = bn - bd - 20;"),
    ("ratio-cut-exponent", W, "dd_core(heap, contents(cw), nn, den, nd, e + sh);", "dd_core(heap, contents(cw), nn, den, nd, e);"),
    ("shift-top-limb", W, "if q + i + 1 < na {", "if q + i + 1 < na - 1 {"),
    ("shift-carry", W, "v = v | ((a[q + i + 1] << (26 - r)) & 67108863);", "v = v | ((a[q + i + 1] << (25 - r)) & 67108863);"),
    ("tail-weights-sum", W, "let (ah, al) = dd.add(sw_h, sw_l, wh[i], wl[i]);", "let (ah, al) = dd.add(sw_h, sw_l, wh[0], wl[0]);"),
    ("tail-divide-by-w", W, "let (b2h, b2l) = dd.div(bh, bl, sw_h, sw_l);", "let (b2h, b2l) = dd.div(bh, bl, 1.0, 0.0);"),
    ("tail-t-complement", W, "dd.add(1.0, 0.0, 0.0 - uh, 0.0 - ul);", "dd.add(0.0, 0.0, 0.0 - uh, 0.0 - ul);"),
    ("tail-t-divisor", W, "float_of(nm1[i]), 0.0);", "float_of(nm1[i] + 1), 0.0);"),
    ("tail-coefficient", W, "dd.div(float_of(2 * (k - 2)), 0.0, kk, 0.0);", "dd.div(float_of(2 * (k - 1)), 0.0, kk, 0.0);"),
    ("tail-denominator-one", W, "dd.add(1.0, 0.0, xh, xl);", "dd.add(0.0, 0.0, xh, xl);"),
    ("tail-between-df", W, "dd.div(b2h, b2l, float_of(k - 1), 0.0);", "dd.div(b2h, b2l, float_of(k), 0.0);"),
    ("tail-df2-three", W, "dd.mul(3.0, 0.0, t_h, t_l);", "dd.mul(2.0, 0.0, t_h, t_l);"),
    ("tail-df2-numerator", W, "dd.div(kk, 0.0, th3, tl3);", "dd.div(float_of(k * k), 0.0, th3, tl3);"),
    ("pair-first-sign", W, "central.put(p, q, t[3 * c..4 * c], l1 + 3, sga < 0);", "central.put(p, q, t[3 * c..4 * c], l1 + 3, sga >= 0);"),
    ("pair-second-sign", W, "central.put(p, q, t[3 * c..4 * c], l1 + 3, sgb >= 0);", "central.put(p, q, t[3 * c..4 * c], l1 + 3, sgb < 0);"),
    ("pair-first-count", W, "welch.mulf(t, 3, 0, l1, f, 9, 3);", "welch.mulf(t, 3, 0, l1, f, 8, 3);"),
    ("pair-second-count", W, "welch.mulf(t, 3, 1, l1, f, 8, 3);", "welch.mulf(t, 3, 1, l1, f, 9, 3);"),
    ("pair-denominator", W, "welch.smul(f, 10, 8, 3, 9, 3);", "welch.smul(f, 10, 8, 3, 8, 3);"),
    ("pair-second-group", W, "psum.magnitude(ab[0..s1], 1, t[c..2 * c])", "psum.magnitude(aa[0..s1], 1, t[c..2 * c])"),
    ("pairs-weights", W, "dd.mul(whv[x], wlv[x], whv[y], wlv[y]);", "dd.mul(whv[x], wlv[x], whv[x], wlv[x]);"),
    ("pairs-square", W, "let (z2h, z2l) = dd.mul(zh, zl, zh, zl);", "let (z2h, z2l) = dd.mul(zh, zl, 1.0, 0.0);"),
    ("pairs-second-from", W, "                                                    var y = x + 1;\n                                                    while y < kk {\n                                                        let gx = giv[x];", "                                                    var y = x;\n                                                    while y < kk {\n                                                        let gx = giv[x];"),
    ("welch-group-zero", W, "if bigint.is_zero(t[2 * c..3 * c], c) {\n                                                wzero = true;", "if bigint.is_zero(t[1 * c..2 * c], c) {\n                                                wzero = true;"),
    ("welch-weight-numerator", W, "welch.smul(f, 6, 0, 3, 5, 6);\n                                                let bwi", "welch.smul(f, 6, 0, 3, 4, 3);\n                                                let bwi"),
    ("welch-counts", W, "nmv[a_] = ns[i] - 1;", "nmv[a_] = ns[i];"),
    ("welch-single", W, "if ns[i] == 1 {\n            single = true;", "if ns[i] == 0 {\n            single = true;"),
    ("welch-two-groups", W, "if kk < 2 || single {", "if kk < 2 {"),
    ("welch-range-weights", W, "if maxw - minw > 360 ||", "if maxw - minw > 500 ||"),
    ("welch-range-distances", W, "maxd - mind > 360 &&", "maxd - mind > 700 &&"),
    ("welch-scale-weights", W, "let ew = 60 - maxw;", "let ew = 50 - maxw;"),
    ("welch-scale-distances", W, "ed = 60 - maxd;", "ed = 50 - maxd;"),
    ("welch-unscale", W, "wf = math.ldexp(fh + fl, 0 - ew - 2 * ed);", "wf = math.ldexp(fh + fl, 0 - ew - ed);"),
    ("welch-overflow", W, "if wf > 1.7e308 {", "if wf > 1.0e30 {"),
    # the small factors and the exact quantities
    ("prep-n-minus-one", A, "central.limbs_of(n - 1, f[4 * welch.slot()..4 * welch.slot() + 3]);", "central.limbs_of(n, f[4 * welch.slot()..4 * welch.slot() + 3]);"),
    ("prep-first-sum", A, "let sg1 = psum.magnitude(a[0..s1], 1, t[0..c]);", "let sg1 = psum.magnitude(a[0..s1], 1, t[c..2 * c]);"),
    ("fold-numerator", A, "welch.mulf(t, 3, 5, welch.len_of(t, 5), f, 0, 3);", "welch.mulf(t, 3, 5, welch.len_of(t, 5), f, 1, 3);"),
    ("fold-product", A, "welch.mulr(t, 4, 2, welch.len_of(t, 2), 6, welch.len_of(t, 6));", "welch.mulr(t, 4, 2, welch.len_of(t, 2), 2, welch.len_of(t, 2));"),
    ("fold-denominator", A, "welch.mulf(t, 3, 6, welch.len_of(t, 6), f, 0, 3);", "welch.mulf(t, 3, 6, welch.len_of(t, 6), f, 1, 3);"),
    ("fold-s2-total", A, "bigint.add_in(t[8 * c..9 * c], c, t[c..c + l2], l2);", "bigint.add_in(t[8 * c..9 * c], c, t[0..l2], l2);"),
    ("total-s1-sign", A, "central.put(pp, qq, t[0..c], l1, sg1 < 0);", "central.put(pp, qq, t[0..c], l1, sg1 >= 0);"),
    ("ssw-times-n", A, "welch.mulf(t, 10, 5, welch.len_of(t, 5), f, 1, 3);", "welch.mulf(t, 10, 5, welch.len_of(t, 5), f, 3, 3);"),
    ("sst-times-d", A, "welch.mulr(t, 11, 9, welch.len_of(t, 9), 6, ld);", "welch.mulr(t, 11, 9, welch.len_of(t, 9), 5, ld);"),
    ("ssb-subtraction", A, "bigint.sub_in(t[12 * c..13 * c], c, t[10 * c..11 * c], c);", "bigint.sub_in(t[12 * c..13 * c], c, t[11 * c..12 * c], c);"),
    ("f-within-df", A, "welch.mulf(t, 3, 12, welch.len_of(t, 12), f, 2, 3);\n                                    welch.mulf(t, 4, 10", "welch.mulf(t, 3, 12, welch.len_of(t, 12), f, 3, 3);\n                                    welch.mulf(t, 4, 10"),
    ("f-between-df", A, "welch.mulf(t, 4, 10, welch.len_of(t, 10), f, 3, 3);", "welch.mulf(t, 4, 10, welch.len_of(t, 10), f, 2, 3);"),
    ("f-variance-zero", A, "let ssw_zero = bigint.is_zero(t[10 * c..11 * c], c);", "let ssw_zero = bigint.is_zero(t[12 * c..13 * c], c);"),
    ("eta-denominator", A, "t[12 * c..13 * c], welch.len_of(t, 12), t[11 * c..12 * c], welch.len_of(t, 11), 0, false);", "t[12 * c..13 * c], welch.len_of(t, 12), t[10 * c..11 * c], welch.len_of(t, 10), 0, false);"),
    ("eta-undefined", A, "if !sst_zero {", "if sst_zero {"),
    ("omega-first-term", A, "central.put(p, q, t[3 * c..4 * c], welch.len_of(t, 3), false);", "central.put(p, q, t[3 * c..4 * c], welch.len_of(t, 3), true);"),
    ("omega-second-term", A, "central.put(p, q, t[3 * c..4 * c], welch.len_of(t, 3), true);\n                                    let so", "central.put(p, q, t[3 * c..4 * c], welch.len_of(t, 3), false);\n                                    let so"),
    ("omega-denominator", A, "welch.sumr(t, 14, 3, 10);", "welch.sumr(t, 14, 3, 12);"),
    ("omega-sign", A, "t[14 * c..15 * c], welch.len_of(t, 14), 0, so < 0);", "t[14 * c..15 * c], welch.len_of(t, 14), 0, so > 0);"),
    # group statistics and Welch's guards
    ("mean-sign", A, "central.ratio_bits(heap, t[0..c], l1, f[0..3], 3, 0 - 1074, sg1 < 0);", "central.ratio_bits(heap, t[0..c], l1, f[0..3], 3, 0 - 1074, sg1 > 0);"),
    ("sd-needs-two", A, "if n >= 2 {\n                                        let (sb, sc)", "if n >= 1 {\n                                        let (sb, sc)"),
    ("k-counts-nonempty", A, "if ns[i] > 0 {\n            k = k + 1;", "if ns[i] > 1 {\n            k = k + 1;"),
    ("too-few-rows", A, "if k < 2 || nn - k < 1 {", "if k < 2 || nn - k < 0 {"),
    # the entries
    ("entry-pooled-within", E, "df2 = float_of(r.nn - r.k);", "df2 = float_of(r.nn - r.k + 1);"),
    ("entry-between-df", E, "float_of(r.k - 1), df2);", "float_of(r.k), df2);"),
    ("entry-pooled-reason", E, "    } else if pooled {\n        if r.status == 2 {", "    } else if !pooled {\n        if r.status == 2 {"),
    ("entry-welch-range", E, "reason = \"welch-range\";", "reason = \"variance-zero\";"),
    ("entry-beyond", E, "} else if wr.status == 4 {", "} else if wr.status == 5 {"),
    ("entry-df-null", E, "if r.status != 1 && (pooled || have) {", "if r.status != 1 {"),
    ("entry-within-pooled", E, "w = json.put_int(heap, w, r.nn - r.k);", "w = json.put_int(heap, w, r.nn);"),
    ("entry-eta-name", E, "w = json.put_string(heap, w, \"eta_squared\");", "w = json.put_string(heap, w, \"eta2\");"),
    ("entry-omega-null", E, "if r.os == 0 && r.oc == 0 {\n        w = json.put_float(heap, w, float_of_bits(r.ob));", "if r.os == 1 && r.oc == 0 {\n        w = json.put_float(heap, w, float_of_bits(r.ob));"),
    ("entry-mean-empty", E, "if counts[i] > 0 {\n            w = json.put_float(heap, w, float_of_bits(gs[3 * i]));", "if counts[i] > 1 {\n            w = json.put_float(heap, w, float_of_bits(gs[3 * i]));"),
    ("entry-sd-one", E, "if counts[i] > 1 {\n            w = json.put_float(heap, w, float_of_bits(gs[3 * i + 1]));", "if counts[i] > 0 {\n            w = json.put_float(heap, w, float_of_bits(gs[3 * i + 1]));"),
    ("entry-sd-index", E, "float_of_bits(gs[3 * i + 1])", "float_of_bits(gs[3 * i])"),
    ("entry-method", E, "    w = json.put_key(heap, w, \"method\");\n    if pooled {", "    w = json.put_key(heap, w, \"method\");\n    if !pooled {"),
    ("entry-p-method", E, "w = json.put_string(heap, w, \"f\");", "w = json.put_string(heap, w, \"t\");"),
    # the driver
    ("main-group-required", F, "\"anova\") && !cli.has(parsed, table, \"group\") {", "\"anova\") && false {"),
    ("main-dispatch", F, "is_test && bytes.equal(cli.operand(args, parsed, 1), \"anova\") {", "is_test && false {"),
    ("run-levels-min", R, "if clean && k < 2 {", "if clean && k < 3 {"),
    ("run-memory-limit", R, "} else if k * nsel > 8192 {", "} else if k * nsel >= 8192 {"),
    ("run-first-group-header", R, "if clean && i == 0 {", "if clean && i == 1 {"),
    ("run-null-rows", R, "var rows_total = gr.nulls;", "var rows_total = 0;"),
    ("run-rows", R, "rows_total = rows_total + info.rows;", "rows_total = rows_total + info.rows + 1;"),
    ("run-numeric-bad", R, "numeric = contents(tn2)[j] > 0 && contents(tb2)[j] == 0;", "numeric = contents(tn2)[j] > 0;"),
    ("run-numeric-total", R, "numeric = contents(tn2)[j] > 0 && contents(tb2)[j] == 0;", "numeric = contents(tb2)[j] == 0;"),
    ("run-accumulate-bad", R, "totbad[j] = totbad[j] + contents(cols.bad)[j];", "totbad[j] = totbad[j] + 0;"),
    ("run-accumulate-numeric", R, "totnum[j] = totnum[j] + contents(cols.numeric)[j];", "totnum[j] = totnum[j] + 0;"),
    ("run-slot", R, "let at = (si * k + i) * g;", "let at = (si * k + i) * g + 1;"),
    ("run-column-slice", R, "if !pooled {\n                                                        wr = welchf", "if pooled {\n                                                        wr = welchf"),
    ("run-level-order", R, "let key = map.key_at(gr.m, contents(gr.order)[i]);", "let key = map.key_at(gr.m, contents(gr.order)[0]);"),
    ("run-header-levels", R, "w = json.put_string(heap, w, map.key_at(gr.m, contents(gr.order)[li]));", "w = json.put_string(heap, w, map.key_at(gr.m, contents(gr.order)[0]));"),
    ("run-adjustment", R, "w = json.put_string(heap, w, cli.text(args, parsed, table, \"correct\"));", "w = json.put_string(heap, w, \"none\");"),
]

EXPLAINED = {
    "welch-weight-numerator": "equivalent: that product only estimates the bit length of a weight to choose the power of two it is scaled by; a few bits of difference move the scale, not the answer (the scaled weights are formed from their own product)",
    "total-s1-sign": "equivalent: the total S1 enters only as S1^2 (LT = N S2 - S1^2); negating every group's sign before summing negates the total",
    "welch-scale-weights": "equivalent: a power of two is exact in double-double, so the scale 2^50 against 2^60 changes the answer only at the edges of the range, and the status 3 (range) is gated at 360 bits, far from them",
    "welch-scale-distances": "equivalent: as welch-scale-weights",
    "pairs-second-from": "equivalent: the pair (x, x) has a zero difference and is skipped by the test for a zero numerator",
    "sd-needs-two": "equivalent in the answer: the entry prints null for a group of fewer than two values whatever this field holds; the guard only keeps the division by n (n - 1) = 0 from being attempted",
    "shift-top-limb": "equivalent: the guard is an index bound; the top limb of a numerator is always zero (the magnitudes have room for 64 bits more than any sum), so what it would contribute is zero",
    "entry-beyond": "unreachable for Welch (argument, not measured): with the weights within 2^360 and the values doubles, F* stays below about 2^466 = 1e140, so the status 4 cannot be set; the pooled F beyond the largest double is reachable and gated",
}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_anova.py"), "--bin", binary, "--files", str(files)], 1200),
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
        print("mutants_anova: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 120):
        print("mutants_anova: the unmutated build fails the gates: refusing to count kills")
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
            dead = not gates(str(d / "build" / "stats"), 120)
        if dead: killed += 1; print("%-32s killed" % n)
        elif n in EXPLAINED: print("%-32s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-32s SURVIVED" % n)
    total = len(M) if only is None else len([1 for n, f, o, w in M if n in only])
    print("mutants_anova: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
