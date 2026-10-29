#!/usr/bin/env python3
"""Mutants of `stats ci boot` (docs/resample.md section 12): lib/stest/boots.cho (the draws, the exact statistics), lib/csvx/boot.cho (the flags, the positions, the entries) and the boot parts of tools/stats/main.cho. Each is a
single edit of the source that must make a gate fail:
    python3 scripts/mutants_boot.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
B = "lib/stest/boots.cho"
D = "lib/csvx/boot.cho"
F = "tools/stats/main.cho"
M = [
    ("key-negative", B, "fn bits_of_key(k: int) -> [] int {\n    if k < 0 {", "fn bits_of_key(k: int) -> [] int {\n    if k < 1 {"),
    ("stat-mean", B, "                        if stat == 0 {\n                            let (mb, mc)", "                        if stat == 1 {\n                            let (mb, mc)"),
    ("sum-squares", B, "                    if stat >= 2 {\n                        psum.add", "                    if stat >= 3 {\n                        psum.add"),
    ("variance-rows", B, "                        } else if n < 2 {\n                            code = 1;", "                        } else if n < 1 {\n                            code = 1;"),
    ("variance-divisor", B, "central.limbs_of(n * (n - 1), contents(dw));", "central.limbs_of(n * n, contents(dw));"),
    ("variance-exponent", B, "0 - 2148, false);\n                                        bits = vb;", "0 - 2147, false);\n                                        bits = vb;"),
    ("mean-exponent", B, "contents(wn), 3, 0 - 1074, neg1);\n                            bits = mb;", "contents(wn), 3, 0 - 1073, neg1);\n                            bits = mb;"),
    ("sd-branch", B, "if stat == 3 && code == 0 {", "if stat == 2 && code == 0 {"),
    ("median-type", B, "quantile.at(heap, contents(kw), n, 7, 5, 1)", "quantile.at(heap, contents(kw), n, 1, 5, 1)"),
    ("median-p", B, "quantile.at(heap, contents(kw), n, 7, 5, 1)", "quantile.at(heap, contents(kw), n, 7, 4, 1)"),
    ("median-branch", B, "    if stat == 1 {\n        return median_bits(heap, picked, n);", "    if stat == 2 {\n        return median_bits(heap, picked, n);"),
    ("draw-kind", B, "rng.derive(rng.derive(seed, 1), b)", "rng.derive(rng.derive(seed, 2), b)"),
    ("draw-replicate", B, "rng.derive(rng.derive(seed, 1), b)", "rng.derive(rng.derive(seed, 1), b + 1)"),
    ("draw-index", B, "wrapping_mul(j + 1, rng.gold())", "wrapping_mul(j, rng.gold())"),
    ("draw-modulus", B, "picked[j] = vals[rng.uniform(z, n)];", "picked[j] = vals[rng.uniform(z, n - 1)];"),
    ("draw-seed", B, "let sb = rng.derive(rng.derive(seed, 1), b);", "let sb = rng.derive(rng.derive(seed + 1, 1), b);"),
    ("replicates-count", B, "while b < reps && code == 0 {", "while b < reps - 1 && code == 0 {"),
    ("replicates-sort", B, "rank.sort_range(out, contents(tw), contents(cw), 0, reps);", "rank.sort_range(out, contents(tw), contents(cw), 0, reps - 1);"),
    ("replicates-store", B, "out[b] = bits_of_key(bits);", "out[b] = bits;"),
    ("seed-length", D, "if len(text) == 0 || len(text) > 19 {", "if len(text) == 0 || len(text) > 18 {"),
    ("seed-empty", D, "if len(text) == 0 || len(text) > 19 {", "if len(text) > 19 {"),
    ("seed-digit-low", D, "if c < 48 || c > 57 {", "if c < 49 || c > 57 {"),
    ("seed-digit-high", D, "if c < 48 || c > 57 {", "if c < 48 || c > 58 {"),
    ("seed-ceiling", D, "if v > (4611686018427387903 - d) / 10 {", "if v > (4611686018427387903 - d) / 9 {"),
    ("seed-required", D, "    if !cli.has(parsed, table, \"seed\") {", "    if cli.has(parsed, table, \"seed\") {"),
    ("replicates-low", D, "if reps < 100 || reps > 100000 {", "if reps < 101 || reps > 100000 {"),
    ("replicates-high", D, "if reps < 100 || reps > 100000 {", "if reps < 100 || reps >= 100000 {"),
    ("stat-median", D, "if bytes.equal(st, \"median\") {\n        stat = 1;", "if bytes.equal(st, \"medianx\") {\n        stat = 1;"),
    ("stat-variance", D, "} else if bytes.equal(st, \"variance\") {\n        stat = 2;", "} else if bytes.equal(st, \"variance\") {\n        stat = 3;"),
    ("stat-sd", D, "} else if bytes.equal(st, \"sd\") {\n        stat = 3;", "} else if bytes.equal(st, \"sd\") {\n        stat = 2;"),
    ("name-median", D, "    if stat == 1 {\n        return \"median\";", "    if stat == 2 {\n        return \"median\";"),
    ("name-sd", D, "    if stat == 3 {\n        return \"sd\";", "    if stat == 2 {\n        return \"sd\";"),
    ("position-low", D, "lo_pos = reps * (scale - lnum) / (2 * scale);", "lo_pos = reps * (scale - lnum) / scale;"),
    ("position-high", D, "hi_pos = (reps * (scale + lnum) + 2 * scale - 1) / (2 * scale) - 1;", "hi_pos = (reps * (scale + lnum)) / (2 * scale) - 1;"),
    ("position-high-off", D, "hi_pos = (reps * (scale + lnum) + 2 * scale - 1) / (2 * scale) - 1;", "hi_pos = (reps * (scale + lnum) + 2 * scale - 1) / (2 * scale);"),
    ("too-few", D, "    if count < 2 {\n        reason = \"too-few-rows\";", "    if count < 3 {\n        reason = \"too-few-rows\";"),
    ("estimate-rows", D, "    if count >= 1 {\n        let (eb, ec)", "    if count >= 2 {\n        let (eb, ec)"),
    ("estimate-code", D, "        if ec == 0 {\n            est = float_of_bits(eb);", "        if ec != 0 {\n            est = float_of_bits(eb);"),
    ("missing", D, "w = json.put_int(heap, w, rows - count);", "w = json.put_int(heap, w, rows);"),
    ("echo-seed", D, "w = json.put_key(heap, w, \"seed\");\n                        w = json.put_int(heap, w, seed);", "w = json.put_key(heap, w, \"seed\");\n                        w = json.put_int(heap, w, reps);"),
    ("echo-replicates", D, "w = json.put_key(heap, w, \"replicates\");\n                        w = json.put_int(heap, w, reps);", "w = json.put_key(heap, w, \"replicates\");\n                        w = json.put_int(heap, w, reps - 1);"),
    ("echo-method", D, "w = json.put_string(heap, w, \"bootstrap-percentile\");", "w = json.put_string(heap, w, \"percentile\");"),
    ("echo-low", D, "lo = float_of_bits(boots.key_of_bits(contents(ow)[lo_pos]));", "lo = float_of_bits(boots.key_of_bits(contents(ow)[lo_pos + 1]));"),
    ("echo-high", D, "hi = float_of_bits(boots.key_of_bits(contents(ow)[hi_pos]));", "hi = float_of_bits(boots.key_of_bits(contents(ow)[hi_pos - 1]));"),
    ("gather-missing", D, "if k != rank.min_int() {\n            out[m] = k;", "if k == rank.min_int() {\n            out[m] = k;"),
    ("selected", D, "                                                if contents(sr2)[j] == 1 {\n                                                    let m = gather", "                                                if contents(sr2)[j] == 0 {\n                                                    let m = gather"),
    ("work-limit", D, "reps * cnt > cli.nat(args, parsed, table, \"max-work\")", "reps * cnt >= cli.nat(args, parsed, table, \"max-work\")"),
    ("work-rows", D, "contents(sr)[j] == 1 && cnt >= 2 && reps * cnt", "contents(sr)[j] == 1 && cnt >= 1 && reps * cnt"),
    ("work-product", D, "reps * cnt > cli.nat(args, parsed, table, \"max-work\")", "reps + cnt > cli.nat(args, parsed, table, \"max-work\")"),
    ("column-pick", D, "pick.select(buffer.bytes(nb0), eb0, nc, cli.text(args, parsed, table, \"cols\"), contents(sw0));", "pick.select(buffer.bytes(nb0), eb0, nc, \"\", contents(sw0));"),
    ("main-kind", F, "&& !bytes.equal(cli.operand(args, parsed, 1), \"corr\") && !bytes.equal(cli.operand(args, parsed, 1), \"boot\") {", "&& !bytes.equal(cli.operand(args, parsed, 1), \"corr\") && !bytes.equal(cli.operand(args, parsed, 1), \"bootx\") {"),
    ("main-cols", F, "\"boot\") && !cli.has(parsed, table, \"cols\") {", "\"boot\") && cli.has(parsed, table, \"cols\") {"),
    ("main-flags", F, "    if is_ci && cli.operand_count(parsed) >= 2 && bytes.equal(cli.operand(args, parsed, 1), \"boot\") {\n        let (e5,", "    if false && is_ci && cli.operand_count(parsed) >= 2 && bytes.equal(cli.operand(args, parsed, 1), \"boot\") {\n        let (e5,"),
    ("main-dispatch", F, "} else if is_ci && bytes.equal(cli.operand(args, parsed, 1), \"boot\") {\n                        let (w2, e2) = boot.boot(", "} else if is_ci && bytes.equal(cli.operand(args, parsed, 1), \"bootx\") {\n                        let (w2, e2) = boot.boot("),
    ("main-stat", F, "boot_seed, boot_reps, boot_stat, lv_num", "boot_seed, boot_reps, 0, lv_num"),
    ("main-level", F, "boot_seed, boot_reps, boot_stat, lv_num, lv_k, float_of(lv_num) / float_of(quantile.pow10(lv_k))", "boot_seed, boot_reps, boot_stat, lv_num, lv_k, 0.95"),
]

EXPLAINED = {}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_ci_boot.py"), "--bin", binary, "--files", str(files)], 1500),
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
        print("mutants_boot: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 40):
        print("mutants_boot: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_boot: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
