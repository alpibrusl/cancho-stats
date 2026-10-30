#!/usr/bin/env python3
"""Mutants of `stats resample perm` (docs/resample.md section 15): lib/stest/permute.cho (the shuffle, the exact comparison), lib/csvx/perm.cho (the flags, the entries) and the resample parts of tools/stats/main.cho. Each is a
single edit of the source that must make a gate fail:
    python3 scripts/mutants_perm.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
PM = "lib/stest/permute.cho"
PE = "lib/csvx/perm.cho"
F = "tools/stats/main.cho"
M = [
    ("key-negative", PM, "fn bits_of_key(k: int) -> [] int {\n    if k < 0 {", "fn bits_of_key(k: int) -> [] int {\n    if k < 1 {"),
    ("dev-sum-sign", PM, "central.put(pos, neg, prod, l1 + 3, ss < 0);", "central.put(pos, neg, prod, l1 + 3, ss > 0);"),
    ("dev-k-sign", PM, "central.put(pos, neg, kmag, l1 + 3, ksign > 0);", "central.put(pos, neg, kmag, l1 + 3, ksign < 0);"),
    ("dev-sum-guard", PM, "    if ss != 0 {\n        bigint.mul(prod", "    if ss > 0 {\n        bigint.mul(prod"),
    ("dev-k-guard", PM, "    if ksign != 0 {\n        central.put", "    if ksign > 0 {\n        central.put"),
    ("m-limbs", PM, "central.limbs_of(m, contents(mw));", "central.limbs_of(m - 1, contents(mw));"),
    ("na-limbs", PM, "central.limbs_of(na, contents(naw));", "central.limbs_of(na + 1, contents(naw));"),
    ("denominator", PM, "central.limbs_of(na * (m - na), contents(dw));", "central.limbs_of(na * m, contents(dw));"),
    ("t-exponent", PM, "0 - 1074, s0 < 0);", "0 - 1073, s0 < 0);"),
    ("t-sign", PM, "0 - 1074, s0 < 0);", "0 - 1074, s0 > 0);"),
    ("total-sum", PM, "sum_first(vals, m, contents(aw)[0..psum.size(1)], contents(ww));", "sum_first(vals, m - 1, contents(aw)[0..psum.size(1)], contents(ww));"),
    ("observed-sum", PM, "sum_first(vals, na, contents(aw)[0..psum.size(1)], contents(ww));", "sum_first(vals, na - 1, contents(aw)[0..psum.size(1)], contents(ww));"),
    ("k-product", PM, "bigint.mul(contents(kw), contents(smw), l1, contents(naw), 3);", "bigint.mul(contents(kw), contents(smw), l1, contents(mw), 3);"),
    ("exact-compare", PM, "                                                                    if bigint.compare(contents(os), central.cap(), contents(o0), central.cap()) >= 0 {\n                                                                        count = count + 1;\n                                                                    }\n                                                                    total = total + 1;", "                                                                    if bigint.compare(contents(os), central.cap(), contents(o0), central.cap()) > 0 {\n                                                                        count = count + 1;\n                                                                    }\n                                                                    total = total + 1;"),
    ("mc-compare", PM, "                                                                    if bigint.compare(contents(os), central.cap(), contents(o0), central.cap()) >= 0 {\n                                                                        count = count + 1;\n                                                                    }\n                                                                    b = b + 1;", "                                                                    if bigint.compare(contents(os), central.cap(), contents(o0), central.cap()) > 0 {\n                                                                        count = count + 1;\n                                                                    }\n                                                                    b = b + 1;"),
    ("shuffle-length", PM, "while d < m - 1 {", "while d < m - 2 {"),
    ("shuffle-top", PM, "let top = m - 1 - d;", "let top = m - d;"),
    ("shuffle-draw", PM, "wrapping_mul(d + 1, rng.gold())", "wrapping_mul(d, rng.gold())"),
    ("shuffle-modulus", PM, "rng.uniform(z, top + 1)", "rng.uniform(z, top)"),
    ("shuffle-kind", PM, "rng.derive(seed, 2)", "rng.derive(seed, 1)"),
    ("shuffle-replicate", PM, "let sb = rng.derive(s1, b);", "let sb = rng.derive(s1, b + 1);"),
    ("shuffle-fresh", PM, "contents(xw)[i] = vals[i];", "contents(xw)[i] = vals[0];"),
    ("shuffle-swap", PM, "contents(xw)[j] = tmp;", "contents(xw)[j] = contents(xw)[top];"),
    ("shuffle-group-a", PM, "sum_first(contents(xw), na, ", "sum_first(contents(xw), m - na, "),
    ("enumeration-end", PM, "while k >= 0 && contents(iw)[k] == m - na + k {", "while k >= 0 && contents(iw)[k] == m - na + k - 1 {"),
    ("enumeration-reset", PM, "contents(iw)[u] = contents(iw)[u - 1] + 1;", "contents(iw)[u] = contents(iw)[u - 1] + 2;"),
    ("mc-total", PM, "total = perms;", "total = perms + 1;"),
    ("flags-low", PE, "if !ok || v < 100 || v > 1000000 {", "if !ok || v < 101 || v > 1000000 {"),
    ("flags-high", PE, "if !ok || v < 100 || v > 1000000 {", "if !ok || v < 100 || v >= 1000000 {"),
    ("flags-digits", PE, "var ok = len(pt) > 0 && len(pt) <= 7;", "var ok = len(pt) > 0 && len(pt) <= 6;"),
    ("flags-digit-high", PE, "if c < 48 || c > 57 {", "if c < 48 || c > 58 {"),
    ("seed-required", PE, "if !bytes.equal(pt, \"all\") && !cli.has(parsed, table, \"seed\") {", "if !cli.has(parsed, table, \"seed\") {"),
    ("seed-required-mc", PE, "if !bytes.equal(pt, \"all\") && !cli.has(parsed, table, \"seed\") {", "if false && !cli.has(parsed, table, \"seed\") {"),
    ("statistic-refusal", PE, "if cli.has(parsed, table, \"statistic\") && !bytes.equal(", "if false && !bytes.equal("),
    ("work-mc", PE, "var work = perms * m;", "var work = perms + m;"),
    ("work-exact", PE, "work = cc * na;", "work = cc;"),
    ("work-limit", PE, "if clean && work > cli.nat(args, parsed, table, \"max-work\") {", "if clean && work >= cli.nat(args, parsed, table, \"max-work\") {"),
    ("exact-too-large", PE, "if cc > 1048576 {", "if cc > 1048577 {"),
    ("p-exact", PE, "p = float_of(count) / float_of(total);", "p = float_of(count) / float_of(total + 1);"),
    ("p-mc", PE, "p = float_of(1 + count) / float_of(total + 1);", "p = float_of(count) / float_of(total + 1);"),
    ("mc-error", PE, "sqrt(p * (1.0 - p) / float_of(perms))", "sqrt(p * (1.0 - p) / float_of(perms + 1))"),
    ("mc-error-exact", PE, "if ok && perms >= 0 {", "if ok {"),
    ("n-a", PE, "w = json.put_key(heap, w, \"a\");\n    w = json.put_int(heap, w, na);", "w = json.put_key(heap, w, \"a\");\n    w = json.put_int(heap, w, nb);"),
    ("gather-b", PE, "let mb = gather(sb, info_b.rows, nc, j, contents(pw), ma);", "let mb = gather(sb, info_b.rows, nc, j, contents(pw), 0);"),
    ("gather-missing", PE, "if k != rank.min_int() {\n            out[m] = k;", "if k == rank.min_int() {\n            out[m] = k;"),
    ("missing", PE, "info_a.rows + info_b.rows - mb, seed, perms);", "info_a.rows + info_b.rows - ma, seed, perms);"),
    ("too-few", PE, "if na < 1 || nb < 1 {\n        reason", "if na < 2 || nb < 1 {\n        reason"),
    ("rows", PE, "w = json.put_int(heap, w, info_a.rows + info_b.rows + gr.nulls);", "w = json.put_int(heap, w, info_a.rows + info_b.rows);"),
    ("seed-echo", PE, "w = json.put_key(heap, w, \"seed\");\n                                                        w = json.put_int(heap, w, seed);", "w = json.put_key(heap, w, \"seed\");\n                                                        w = json.put_int(heap, w, seed + 1);"),
    ("level-a", PE, "w = json.put_string(heap, w, key_a);", "w = json.put_string(heap, w, key_b);"),
    ("exact-name", PE, "w = json.put_string(heap, w, \"exact\");", "w = json.put_string(heap, w, \"monte-carlo\");"),
    ("group-levels", PE, "if clean && nlevels != 2 {", "if clean && nlevels != 3 {"),
    ("main-verb", F, "!bytes.equal(cli.operand(args, parsed, 0), \"resample\") {\n        e = rule.refuse(heap, e, \"stats.bad-verb\"", "!bytes.equal(cli.operand(args, parsed, 0), \"resamplex\") {\n        e = rule.refuse(heap, e, \"stats.bad-verb\""),
    ("main-kind", F, "} else if !bytes.equal(cli.operand(args, parsed, 1), \"perm\") {", "} else if !bytes.equal(cli.operand(args, parsed, 1), \"permx\") {"),
    ("main-operands", F, "} else if cli.operand_count(parsed) != 3 {\n            e = rule.refuse(heap, e, \"args.missing-operand\", \"resample names", "} else if cli.operand_count(parsed) > 3 {\n            e = rule.refuse(heap, e, \"args.missing-operand\", \"resample names"),
    ("main-required", F, "} else if !cli.has(parsed, table, \"cols\") || !cli.has(parsed, table, \"group\") {\n            e = rule.refuse(heap, e, \"stats.test-columns\", \"resample perm", "} else if !cli.has(parsed, table, \"cols\") && !cli.has(parsed, table, \"group\") {\n            e = rule.refuse(heap, e, \"stats.test-columns\", \"resample perm"),
    ("main-dispatch", F, "} else if is_resample {\n                        let (w2, e2)", "} else if false {\n                        let (w2, e2)"),
    ("main-flags", F, "if is_resample && cli.operand_count(parsed) >= 2 && bytes.equal(cli.operand(args, parsed, 1), \"perm\") {", "if false && is_resample && cli.operand_count(parsed) >= 2 && bytes.equal(cli.operand(args, parsed, 1), \"perm\") {"),
    ("main-seed", F, "perm_seed, perm_n, w, e);", "0, perm_n, w, e);"),
    ("main-permutations", F, "perm_seed, perm_n, w, e);", "perm_seed, 100, w, e);"),
]

EXPLAINED = {}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_resample_perm.py"), "--bin", binary, "--files", str(files)], 1500),
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
        print("mutants_perm: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 40):
        print("mutants_perm: the unmutated build fails the gates: refusing to count kills")
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
            dead = not gates(str(d / "build" / "stats"), 40)
        if dead: killed += 1; print("%-32s killed" % n)
        elif n in EXPLAINED: print("%-32s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-32s SURVIVED" % n)
    total = len(M) if only is None else len([1 for n, f, o, w in M if n in only])
    print("mutants_perm: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
