#!/usr/bin/env python3
"""Mutants of `stats test fisher` (docs/rank.md R-H): lib/stest/fisher.cho (the weights, the sums, the odds ratio), lib/stest/answer10.cho (the entries), the fisher parts of tools/stats/main.cho. Each is applied to a copy, the program rebuilt, scripts/gate_fisher.py and the refusals gate rerun; killed = a gate fails or the program traps or hangs. A mutant that does not compile is INVALID. The gates first run on the UNMUTATED build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.
    python3 scripts/mutants_fisher.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
W = "lib/stest/fisher.cho"
E = "lib/stest/answer10.cho"
R = "lib/csvx/run.cho"
F = "tools/stats/main.cho"
K = "lib/csvx/cat.cho"
M = [
    ("limit", W, "    return 10000;", "    return 10001;"),
    ("mul-small", W, "let t = v[i] * m + carry;", "let t = v[i] * m;"),
    ("div-small-base", W, "let t = rem * 67108864 + v[i];", "let t = rem * 67108863 + v[i];"),
    ("binomial-step", W, "len = mul_small(v, len, n - kk + i);", "len = mul_small(v, len, n - kk + i + 1);"),
    ("binomial-symmetry", W, "if n - k < kk {", "if n - k > kk {"),
    ("first-weight-branch", W, "if c1 <= r2 {", "if c1 < r2 {"),
    ("step-mul-row", W, "var l = mul_small(w, len, r1 - x);", "var l = mul_small(w, len, r1 - x - 1);"),
    ("step-div-row", W, "l = div_small(w, l, x + 1);", "l = div_small(w, l, x + 2);"),
    ("step-mul-col", W, "l = mul_small(w, l, c1 - x);", "l = mul_small(w, l, c1 - x - 1);"),
    ("step-div-col", W, "l = div_small(w, l, r2 - c1 + x + 1);", "l = div_small(w, l, r2 - c1 + x);"),
    ("xmin", W, "if c1 - r2 > 0 {", "if c1 - r2 > 1 {"),
    ("xmax", W, "if r1 < xmax {", "if r1 > xmax {"),
    ("observed", W, "if x == a {", "if x == b {"),
    ("not-larger", W, "<= 0 {\n                    bigint.add_in(s[2 * cap..3 * cap]", "< 0 {\n                    bigint.add_in(s[2 * cap..3 * cap]"),
    ("sum-all", W, "bigint.add_in(s[3 * cap..4 * cap], cap, s[0..cap], len);", "bigint.add_in(s[3 * cap..4 * cap], cap, s[cap..2 * cap], lobs);"),
    ("sum-not-larger", W, "bigint.add_in(s[2 * cap..3 * cap], cap, s[0..cap], len);", "bigint.add_in(s[2 * cap..3 * cap], cap, s[cap..2 * cap], lobs);"),
    ("second-pass-start", W, "bigint.clear(s[0..cap], cap);", "bigint.clear(s[2 * cap..3 * cap], cap);"),
    ("odds-ratio", W, "orb = bits_of(float_of(num) / float_of(den));", "orb = bits_of(float_of(den) / float_of(num));"),
    ("odds-infinite", W, "if num == 0 {\n                orkind = 1;", "if num == 1 {\n                orkind = 1;"),
    ("empty", W, "if n == 0 {\n        status = 1;", "if n == 1 {\n        status = 1;"),
    ("entry-status", E, "if r.status == 1 {\n        reason = \"too-few-rows\";", "if r.status == 2 {\n        reason = \"too-few-rows\";"),
    ("entry-pc", E, "if r.pc == 0 {\n            have_p = true;", "if r.pc == 1 {\n            have_p = true;"),
    ("entry-odds-defined", E, "if r.orkind == 0 {\n            have_or = true;", "if r.orkind == 1 {\n            have_or = true;"),
    ("entry-reason-swap", E, "if r.orkind == 1 {\n                reason = \"odds-ratio-undefined\";", "if r.orkind == 2 {\n                reason = \"odds-ratio-undefined\";"),
    ("entry-statistic", E, "w = json.put_float(heap, w, float_of_bits(r.orb));\n    } else {\n        w = json.put_null(heap, w);\n    }\n    w = json.put_key(heap, w, \"df\");", "w = json.put_float(heap, w, float_of_bits(r.pb));\n    } else {\n        w = json.put_null(heap, w);\n    }\n    w = json.put_key(heap, w, \"df\");"),
    ("entry-p-method", E, "\"exact, probabilities not larger than the observed\"", "\"exact\""),
    ("entry-effect-name", E, "\"odds_ratio\"", "\"odds\""),
    ("entry-n", E, "w = json.put_int(heap, w, r.n);", "w = json.put_int(heap, w, r.a);"),
    ("entry-missing", E, "w = json.put_int(heap, w, rows - r.n);", "w = json.put_int(heap, w, rows);"),
    ("entry-method", E, "w = json.put_string(heap, w, \"fisher\");", "w = json.put_string(heap, w, \"fisherx\");"),
    ("entry-level-a", E, "map.key_at(d.cross.la, contents(d.order_a)[0])", "map.key_at(d.cross.la, contents(d.order_a)[1])"),
    ("entry-level-vs", E, "map.key_at(d.cross.lb, contents(d.order_b)[1]));\n        w = json.end_array", "map.key_at(d.cross.lb, contents(d.order_b)[0]));\n        w = json.end_array"),
    ("entry-table", E, "w = json.put_int(heap, w, r.b);\n        w = json.put_int(heap, w, r.c);", "w = json.put_int(heap, w, r.c);\n        w = json.put_int(heap, w, r.b);"),
    ("entry-assumptions", E, "if r.status != 1 {\n        w = json.put_key(heap, w, \"levels\");", "if r.status != 2 {\n        w = json.put_key(heap, w, \"levels\");"),
    ("main-levels", K, "} else if dr.rows_ != 2 || dr.cols_ != 2 {", "} else if dr.rows_ != 2 && dr.cols_ != 2 {"),
    ("main-empty", K, "if dr.rows_ == 0 {\n                                w = answer10.put_fisher", "if dr.rows_ == 1 {\n                                w = answer10.put_fisher"),
    ("main-limit", K, "} else if dr.cross.used > fisher.limit() {", "} else if dr.cross.used >= fisher.limit() {"),
    ("main-counts", K, "fisher.run(heap, cn[0], cn[1], cn[2], cn[3])", "fisher.run(heap, cn[0], cn[2], cn[1], cn[3])"),
    ("main-vs-required", F, "\"fisher\") && !cli.has(parsed, table, \"vs\") {", "\"fisher\") && false {"),
    ("main-dispatch", F, "\"fisher\") {\n                        let (w2, e2) = cat.fisher_one(", "\"fisherx\") {\n                        let (w2, e2) = cat.fisher_one("),
    ("main-kind-accepted", F, "&& !bytes.equal(cli.operand(args, parsed, 1), \"fisher\") {", "&& false {"),
    ("main-header", K, "w = json.put_string(heap, w, \"fisher\");\n        w = json.put_key(heap, w, \"vs\");", "w = json.put_string(heap, w, \"chisq\");\n        w = json.put_key(heap, w, \"vs\");"),
]

EXPLAINED = {
    "binomial-symmetry": "equivalent in the answer: C(n, k) = C(n, n - k), so the loop over k instead of min(k, n - k) gives the same integer; only the number of steps changes",
    "first-weight-branch": "equivalent: at c1 = r2 the two branches are C(r2, c1) = 1 and C(r1, 0) = 1",
    "xmax": "equivalent: beyond min(r1, c1) the factor r1 - x or c1 - x of the step is 0, so every later weight is 0 and adds nothing to either sum",
    "second-pass-start": "equivalent: after the first pass the weight slot holds w(xmax), but the binomial writes limb 0 and the multiplication writes a new limb only where it extends the length, and every later operation reads only the limbs in use, so the limbs above are never read",
}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_fisher.py"), "--bin", binary, "--files", str(files)], 1500),
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
        print("mutants_fisher: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 100):
        print("mutants_fisher: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_fisher: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
