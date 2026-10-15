#!/usr/bin/env python3
"""Mutants of tools/stats/main.cho (docs/describe.md S1-G): each is applied to a copy, the driver rebuilt, scripts/gate_describe.py rerun; killed = the gate fails or the driver traps. A mutant that does
not compile is INVALID. Exit 1 if one survives or is invalid.   python3 scripts/mutants_stats.py [--check]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
F = "tools/stats/main.cho"
M = [
    ("mean-sign", "ratio_bits(heap, contents(w1), l1, contents(wn), 3, 0 - 1074, neg1)", "ratio_bits(heap, contents(w1), l1, contents(wn), 3, 0 - 1074, false)"),
    ("mean-unit", "contents(wn), 3, 0 - 1074, neg1)", "contents(wn), 3, 0 - 1073, neg1)"),
    ("variance-unit", "contents(dw), 3, 0 - 2148, false)", "contents(dw), 3, 0 - 2147, false)"),
    ("variance-ddof-ignored", "limbs_of(n * (n - ddof), contents(dw));", "limbs_of(n * n, contents(dw));"),
    ("numerator-sum-not-difference", "bigint.sub_in(contents(aw), l2 + 5, contents(bw), 2 * l1 + 2);", "bigint.add_in(contents(aw), l2 + 5, contents(bw), 2 * l1 + 2);"),
    ("count-limb-shift", "    out[1] = v >> 26 & 67108863;", "    out[1] = v >> 25 & 67108863;"),
    ("min-comparison", "if numeric[j] == 0 || x < mn[j] {", "if numeric[j] == 0 || x > mn[j] {"),
    ("max-comparison", "if numeric[j] == 0 || x > mx[j] {", "if numeric[j] == 0 || x < mx[j] {"),
    ("min-first-value", "if numeric[j] == 0 || x < mn[j] {", "if x < mn[j] {"),
    ("count-twice", "nonmiss[j] = nonmiss[j] + 1;", "nonmiss[j] = nonmiss[j] + 2;"),
    ("missing-from-numeric", "w = json.put_int(heap, w, rows - nonmiss);", "w = json.put_int(heap, w, rows - numeric);"),
    ("quoted-cell-is-a-number", "            if quoted == 1 {\n                bad[j] = 1;\n            } else if bad[j] == 0 {", "            if quoted == 2 {\n                bad[j] = 1;\n            } else if bad[j] == 0 {"),
    ("sd-is-variance", "        sd = sqrt(variance);", "        sd = variance;"),
    ("variance-guard", "if n - ddof >= 1 {", "if n - ddof >= 0 {"),
    ("second-power-accumulator", "psum.add(accs[base + psum.size(1)..base + acc_stride()], 2, bits_of(x), work);", "psum.add(accs[base + psum.size(1)..base + acc_stride()], 3, bits_of(x), work);"),
    ("blank-line-one-column", "} else if len(line) > 0 || info.ncols == 1 {", "} else if len(line) > 0 {"),
    ("ragged-check", "    if n != ncols {\n        return 1;", "    if n < ncols {\n        return 1;"),
    ("row-limit-edge", "if info.rows >= max_rows {", "if info.rows > max_rows {"),
    ("column-limit-edge", "} else if seps + 1 > max_cols {", "} else if seps + 1 >= max_cols {"),
    ("empty-kind", "        if nonmiss == 0 {\n            w = json.put_string(heap, w, \"empty\");", "        if nonmiss == 1 {\n            w = json.put_string(heap, w, \"empty\");"),
    ("verb-name", "!bytes.equal(cli.operand(args, parsed, 0), \"describe\")", "!bytes.equal(cli.operand(args, parsed, 0), \"describe \")"),
    ("delimiter-fixed", "let delim = int_of(cli.text(args, parsed, table, \"delim\")[0]);\n    let max_rows", "let delim = 44;\n    let max_rows"),
    ("ddof-ceiling", "if cli.nat(args, parsed, table, \"ddof\") > 1 {", "if cli.nat(args, parsed, table, \"ddof\") > 2 {"),
    ("float-key-of-negative", "        key = bits ^ 0x7fffffffffffffff;", "        key = bits;"),
    ("multiline-record-taken", "    if inside {\n        return 2;\n    }", ""),
    ("overflow-code-ignored", "                ok = vc == 0;", "                ok = true;"),
]

EXPLAINED = {}


def main():
    t = (ROOT / F).read_text()
    for n, o, w in M:
        if t.count(o) != 1 or o == w:
            print("mutant %s: pattern occurs %d times" % (n, t.count(o)))
            return 1
    if "--check" in sys.argv:
        print("mutants_stats: %d mutants, every pattern applies once" % len(M)); return 0
    bad, killed = [], 0
    for n, o, w in M:
        with tempfile.TemporaryDirectory() as d:
            d = pathlib.Path(d)
            for x in ("lib", "tools", "generated"): shutil.copytree(ROOT / x, d / x)
            shutil.copy(ROOT / "cancho.toml", d / "cancho.toml")
            (d / "build").mkdir(); shutil.copytree(ROOT / "build" / "deps", d / "build" / "deps")
            f = d / F; f.write_text(f.read_text().replace(o, w, 1))
            r = subprocess.run([CANCHO, "build", "--bin", "stats"], cwd=d, capture_output=True, text=True)
            if r.returncode != 0:
                bad.append(n); print("%-28s INVALID: %s" % (n, r.stderr[-120:].replace("\n", " "))); continue
            try:
                g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe.py"), "--driver", str(d / "build" / "stats"), "--files", "60"], capture_output=True, text=True, timeout=300)
                dead = g.returncode != 0
                if not dead:
                    g2 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe_refusals.py"), "--bin", str(d / "build" / "stats"), "--fuzz", "150"], capture_output=True, text=True, timeout=600)
                    dead = g2.returncode != 0
            except subprocess.TimeoutExpired:
                dead = True
        if dead: killed += 1; print("%-28s killed" % n)
        elif n in EXPLAINED: print("%-28s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-28s SURVIVED" % n)
    print("mutants_stats: %d killed of %d, %d survived and explained, %d unexplained" % (killed, len(M), len(M) - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
