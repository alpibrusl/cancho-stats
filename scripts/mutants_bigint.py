#!/usr/bin/env python3
"""Mutants of lib/acc/bigint.cho (docs/describe.md S1-G): each is applied to a copy, the driver rebuilt, scripts/gate_bigint.py rerun; killed = the gate fails or the driver traps. A mutant that does
not compile is INVALID. Exit 1 if one survives or is invalid.   python3 scripts/mutants_bigint.py [--check]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
F = "lib/acc/bigint.cho"
M = [
    ("compare-greater-sign", "        if x > y {\n            return 1;", "        if x < y {\n            return 1;"),
    ("compare-equal-case", "        if x < y {\n            return -1;\n        }\n        i = i - 1;", "        i = i - 1;"),
    ("bit-length-off-by-one", "            return i * 26 + b;", "            return i * 26 + b + 1;"),
    ("bit-at-limb-index", "    return a[pos / 26] >> pos % 26 & 1;", "    return a[pos / 26] >> pos % 25 & 1;"),
    ("set-bit-or-to-xor", "a[pos / 26] = a[pos / 26] | 1 << pos % 26;", "a[pos / 26] = a[pos / 26] ^ 1 << pos % 26;"),
    ("shl-add-carry", "        let t = (a[i] << k) + carry;\n        a[i] = t & mask();\n        carry = t >> 26;", "        let t = (a[i] << k) + carry;\n        a[i] = t & mask();\n        carry = t >> 25;"),
    ("sub-borrow-base", "            t = t + 67108864;", "            t = t + 67108863;"),
    ("sub-borrow-flag", "            owe = 1;\n        } else {", "            owe = 0;\n        } else {"),
    ("add-carry-dropped", "        a[i] = t & mask();\n        carry = t >> 26;\n        i = i + 1;\n    }\n    return carry;\n}\n\n// out[0 .. na + nb)", "        a[i] = t & mask();\n        carry = 0;\n        i = i + 1;\n    }\n    return carry;\n}\n\n// out[0 .. na + nb)"),
    ("mul-carry-row", "        out[i + nb] = out[i + nb] + carry;", "        out[i + nb] = out[i + nb];"),
    ("mul-partial-product", "let t = out[i + j] + a[i] * b[j] + carry;", "let t = out[i + j] + a[i] * b[j];"),
    ("div-shift-window", "        if pos >= shift {\n            b = bit_at(num, nn, pos - shift);", "        if pos > shift {\n            b = bit_at(num, nn, pos - shift);"),
    ("div-compare-strict", "        if compare(r, nd + 1, den, nd) >= 0 {", "        if compare(r, nd + 1, den, nd) > 0 {"),
    ("div-sticky-dropped", "    if is_zero(r, nd + 1) {\n        return 0;\n    }\n    return 1;", "    return 0;"),
    ("round-guard-tie-to-even", "    if guard == 1 && (rest || m & 1 == 1) {", "    if guard == 1 && rest {"),
    ("round-guard-always-up", "    if guard == 1 && (rest || m & 1 == 1) {", "    if guard == 1 {"),
    ("round-subnormal-clamp", "        lsb = -1074 - e;", "        lsb = -1073 - e;"),
    ("round-renormalise", "    if m == 9007199254740992 {\n        m = 4503599627370496;\n        q2 = q2 + 1;\n    }", ""),
    ("round-overflow-edge", "        if biased >= 2047 {", "        if biased > 2047 {"),
    ("round-sticky-ignored", "    var rest = sticky == 1;", "    var rest = false;"),
    ("round-exponent-bias", "        let biased = q2 + 1075;", "        let biased = q2 + 1074;"),
]

EXPLAINED = {"set-bit-or-to-xor": "equivalent: set_bit is called once per quotient bit on a zeroed quotient, and a bit set once is the same under or and xor"}


def main():
    t = (ROOT / F).read_text()
    for n, o, w in M:
        if t.count(o) != 1 or o == w:
            print("mutant %s: pattern occurs %d times" % (n, t.count(o)))
            return 1
    if "--check" in sys.argv:
        print("mutants_bigint: %d mutants, every pattern applies once" % len(M)); return 0
    bad, killed = [], 0
    for n, o, w in M:
        with tempfile.TemporaryDirectory() as d:
            d = pathlib.Path(d)
            for x in ("lib", "tools"): shutil.copytree(ROOT / x, d / x)
            shutil.copy(ROOT / "cancho.toml", d / "cancho.toml")
            (d / "build").mkdir(); shutil.copytree(ROOT / "build" / "deps", d / "build" / "deps")
            f = d / F; f.write_text(f.read_text().replace(o, w, 1))
            r = subprocess.run([CANCHO, "build", "--bin", "bigint-check"], cwd=d, capture_output=True, text=True)
            if r.returncode != 0:
                bad.append(n); print("%-28s INVALID: %s" % (n, r.stderr[-120:].replace("\n", " "))); continue
            try:
                g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_bigint.py"), "--driver", str(d / "build" / "bigint-check"), "--n", "600"], capture_output=True, text=True, timeout=300)
                dead = g.returncode != 0
            except subprocess.TimeoutExpired:
                dead = True
        if dead: killed += 1; print("%-28s killed" % n)
        elif n in EXPLAINED: print("%-28s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-28s SURVIVED" % n)
    print("mutants_bigint: %d killed of %d, %d survived and explained, %d unexplained" % (killed, len(M), len(M) - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
