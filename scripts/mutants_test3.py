#!/usr/bin/env python3
"""Mutants of `stats test chisq` (docs/test.md T-J): lib/csvx/cross.cho (the counting pass and the table), lib/stest/chisq.cho (the statistic, the proportions of --expected), lib/stest/answer3.cho and the chisq part of
tools/stats/main.cho. Each is applied to a copy, the program rebuilt, scripts/gate_test3.py and the refusals gate rerun; killed = a gate fails or the program traps or hangs. A mutant that does not compile is INVALID. The
gates first run on the UNMUTATED build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.   python3 scripts/mutants_test3.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
X = "lib/csvx/cross.cho"
Q = "lib/stest/chisq.cho"
B = "lib/stest/answer3.cho"
F = "tools/stats/main.cho"
M = [
    # the counting pass
    ("second-cell-missing-taken", X, "if tb == fb && cs[3 * ib + 2] == 0 {\n                                        present = false;", "if tb == fb && cs[3 * ib + 2] == 0 {\n                                        present = true;"),
    ("first-cell-missing-taken", X, "var present = !(ta == fa && cs[3 * ia + 2] == 0);", "var present = true;"),
    ("quoted-empty-is-missing", X, "var present = !(ta == fa && cs[3 * ia + 2] == 0);", "var present = !(ta == fa);"),
    ("count-twice", X, "m = map.put(heap, m, buffer.bytes(kr), count + 1);", "m = map.put(heap, m, buffer.bytes(kr), count + 2);"),
    ("first-level-limit-edge", X, "room_a = map.size(lar) < max_groups;", "room_a = map.size(lar) <= max_groups;"),
    ("second-level-limit-edge", X, "room_b = map.size(lbr) < max_groups;", "room_b = map.size(lbr) <= max_groups;"),
    ("second-ambiguity", X, "} else if mb > 1 {", "} else if mb > 2 {"),
    ("second-no-match", X, "if mb == 0 {", "if mb == 1 {"),
    ("first-ambiguity", X, "} else if ma > 1 {", "} else if ma > 2 {"),
    ("second-column-ignored", X, "if two && !stop {", "if false && !stop {"),
    ("pair-key-order", X, "kb = put_key(heap, kb, line[fa..ta], line[fb..tb]);", "kb = put_key(heap, kb, line[fb..tb], line[fa..ta]);"),
    ("table-transposed", X, "cn[ra * c + rb] = map.value_at(mr, en);", "cn[rb * r + ra] = map.value_at(mr, en);"),
    ("levels-descending", X, "map.key_at(lr1, en)) > 0 {", "map.key_at(lr1, en)) < 0 {"),
    ("ranks-all-zero", X, "map.set_value_at(lw, en, k);", "map.set_value_at(lw, en, 0);"),
    ("one-column-width", X, "if !two {\n        c = 1;", "if !two {\n        c = 2;"),
    # the statistic
    ("independence-sign", Q, "n * counts[i * c + j] - rc, n, rc, 0);", "n * counts[i * c + j] + rc, n, rc, 0);"),
    ("expected-sum-not-product", Q, "let rc = rsum[i] * csum[j];", "let rc = rsum[i] + csum[j];"),
    ("minimum-expected", Q, "if first || rc < rc_min {", "if first || rc > rc_min {"),
    ("below-five-edge", Q, "if rc < 5 * n {", "if rc <= 5 * n {"),
    ("below-five-threshold", Q, "if rc < 5 * n {", "if rc < 4 * n {"),
    ("df-independence", Q, "df: (r - 1) * (c - 1)", "df: r * c - 1"),
    ("effect-max-not-min", Q, "    if c < r {\n        m = c;", "    if c > r {\n        m = c;"),
    ("effect-factor", Q, "effect = sqrt(sh / (float_of(n) * float_of(m - 1)));", "effect = sqrt(sh / float_of(n));"),
    ("fit-weights", Q, "let np = n * props[i];", "let np = n * props[0];"),
    ("fit-observed-scale", Q, "counts[i] * d - np, np, d, 0);", "counts[i] - np, np, d, 0);"),
    ("fit-below-five-edge", Q, "if np < 5 * d {", "if np <= 5 * d {"),
    ("fit-effect", Q, "effect = sqrt(sh / float_of(n));\n    }\n    return ChiSq { n: n, rows_: k", "effect = sqrt(sh / float_of(k));\n    }\n    return ChiSq { n: n, rows_: k"),
    ("fit-df", Q, "df: k - 1,", "df: k,"),
    ("int-split-low", Q, "float_of(v >> 31) * 2147483648.0", "float_of(v >> 30) * 2147483648.0"),
    ("term-zero-skipped-wrongly", Q, "if num == 0 {\n        return (sh, sl);", "if num == 1 {\n        return (sh, sl);"),
    ("ratio-empty", Q, "if b == 0 {\n        return 0.0;", "if b == 7 {\n        return 0.0;"),
    # --expected
    ("expected-scale", Q, "scale = scale * 10;", "scale = scale * 9;"),
    ("expected-sum-above-only", Q, "if status == 0 && sum != 1000000000 {", "if status == 0 && sum > 1000000000 {"),
    ("expected-two-points", Q, "} else if ch == 46 && !seen_point {", "} else if ch == 46 {"),
    ("expected-decimals", Q, "if decimals > 9 {", "if decimals > 10 {"),
    ("expected-zero", Q, "if v == 0 {\n                    status = 2;", "if v == 1 {\n                    status = 2;"),
    ("expected-digits-needed", Q, "if digits == 0 {\n                status = 1;", "if digits == 7 {\n                status = 1;"),
    # the answer
    ("v-name", B, "w = json.put_string(heap, w, \"cramers_v\");", "w = json.put_string(heap, w, \"cramer_v\");"),
    ("w-name", B, "w = json.put_string(heap, w, \"cohen_w\");", "w = json.put_string(heap, w, \"cohen_d\");"),
    ("method-names", B, "w = json.put_string(heap, w, \"goodness of fit\");", "w = json.put_string(heap, w, \"fit\");"),
    ("p-method", B, "w = json.put_string(heap, w, \"chi2\");", "w = json.put_string(heap, w, \"t\");"),
    ("rounding-label", B, "w = json.put_string(heap, w, \"double-double\");", "w = json.put_string(heap, w, \"exact\");"),
    ("missing-count", B, "w = json.put_int(heap, w, rows - r.n);", "w = json.put_int(heap, w, rows);"),
    ("levels-first-only", B, "map.key_at(d.cross.la, contents(d.order_a)[k])", "map.key_at(d.cross.la, contents(d.order_a)[0])"),
    ("one-level-not-null", B, "} else if r.df < 1 {\n        reason = \"too-few-levels\";", "} else if r.df < 0 {\n        reason = \"too-few-levels\";"),
    ("empty-reason", B, "reason = \"too-few-rows\";", "reason = \"too-few-levels\";"),
    ("df-printed", B, "if r.n > 0 && r.df >= 1 {", "if r.n > 0 && r.df >= 0 {"),
    ("below-five-printed", B, "w = json.put_int(heap, w, r.below5);", "w = json.put_int(heap, w, r.rows_);"),
    ("table-dims", B, "w = json.put_int(heap, w, r.rows_);\n            w = json.put_key(heap, w, \"cols\");", "w = json.put_int(heap, w, r.cols_);\n            w = json.put_key(heap, w, \"cols\");"),
    # the verb
    ("vs-expected-exclusive", F, "cli.has(parsed, table, \"vs\") == cli.has(parsed, table, \"expected\")", "cli.has(parsed, table, \"vs\") != cli.has(parsed, table, \"expected\")"),
    ("fit-flag", F, "let fit = cli.has(parsed, table, \"expected\");", "let fit = false;"),
    ("level-count-check", F, "if dr.rows_ != kprops {", "if dr.rows_ > kprops {"),
    ("chisq-dispatch", F, "if is_test && bytes.equal(cli.operand(args, parsed, 1), \"chisq\") {", "if false {"),
    ("empty-item", F, "if len(bytes.field(list, ',', i + 1)) == 0 {", "if len(bytes.field(list, ',', i + 1)) == 1 {"),
    ("expected-echo", F, "float_of(contents(pr0)[k]) / 1000000000.0", "float_of(contents(pr0)[k]) / 100000000.0"),
    ("columns-after-first", F, "while ok && i < ncols_named {\n        let name = bytes.field(list, ',', i + 1);\n        match place.open_operand(fs, root, rel, full) {\n            Opened::Failed(reason) => {\n                e = fail.io_error(heap, e, reason, false, rel);\n            }\n            Opened::Ok(opened) => {\n                var file = opened;\n                borrow mut file as &!handle in {\n                    let (x, e2) = cross.cross(heap, handle, rel, name, cli.text(args, parsed, table, \"vs\"), !fit", "while ok && i < 1 {\n        let name = bytes.field(list, ',', i + 1);\n        match place.open_operand(fs, root, rel, full) {\n            Opened::Failed(reason) => {\n                e = fail.io_error(heap, e, reason, false, rel);\n            }\n            Opened::Ok(opened) => {\n                var file = opened;\n                borrow mut file as &!handle in {\n                    let (x, e2) = cross.cross(heap, handle, rel, name, cli.text(args, parsed, table, \"vs\"), !fit"),
    ("kind-accepted", F, "&& !bytes.equal(cli.operand(args, parsed, 1), \"chisq\") && !bytes", "&& false && !bytes"),
]

EXPLAINED = {}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_test3.py"), "--bin", binary, "--files", str(files)], 400),
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
        print("mutants_test3: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 100):
        print("mutants_test3: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_test3: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
