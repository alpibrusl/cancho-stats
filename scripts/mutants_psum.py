#!/usr/bin/env python3
"""Mutants of lib/acc/psum.cho (docs/describe.md S1-G): each is applied to a copy, the driver rebuilt, scripts/gate_psum.py rerun; killed = the gate fails or the driver traps. A mutant that does
not compile is INVALID. Exit 1 if one survives or is invalid.   python3 scripts/mutants_psum.py [--check]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
F = "lib/acc/psum.cho"
M = [
    ("hidden-bit-dropped", "        m = m | 4503599627370496;\n        f = ef - 1;", "        f = ef - 1;"),
    ("exponent-offset", "        f = ef - 1;", "        f = ef;"),
    ("shift-uses-f-not-kf", "return add_shifted(acc, base, len, work, poff, npw, k * f);", "return add_shifted(acc, base, len, work, poff, npw, f);"),
    ("sign-of-even-power", "if bits < 0 && k % 2 == 1 {", "if bits < 0 {"),
    ("sign-of-odd-power", "if bits < 0 && k % 2 == 1 {", "if bits < 0 && k % 2 == 0 {"),
    ("limbs-too-few", "return (k * 2151 + 64 + 25) / 26;", "return (k * 2098 + 25) / 26;"),
    ("mantissa-limb-1", "    work[1] = m >> 26 & mask();", "    work[1] = m >> 27 & mask();"),
    ("mantissa-limb-2", "    work[2] = m >> 52;", "    work[2] = m >> 51;"),
    ("square-operands", "bigint.mul(work[3..9], work[0..3], 3, work[0..3], 3);", "bigint.mul(work[3..9], work[0..3], 3, work[0..3], 2);"),
    ("cube-uses-m-twice", "bigint.mul(work[9..18], work[3..9], 6, work[0..3], 3);", "bigint.mul(work[9..18], work[3..9], 6, work[3..9], 3);"),
    ("fourth-power-operands", "bigint.mul(work[18..30], work[3..9], 6, work[3..9], 6);", "bigint.mul(work[18..30], work[3..9], 6, work[0..3], 3);"),
    ("fourth-power-selected", "    if k == 4 {\n        bigint.mul", "    if k == 5 {\n        bigint.mul"),
    ("carry-high-part", "        carry = (t >> 26) + (v >> 26);", "        carry = (t >> 26);"),
    ("carry-propagation", "    while carry != 0 && j < len {", "    while carry != 0 && j < q {"),
    ("magnitude-compare", "    if c > 0 {\n        bigint.copy", "    if c < 0 {\n        bigint.copy"),
    ("magnitude-sign", "        bigint.sub_in(out, len, acc[len..2 * len], len);\n        return 1;", "        bigint.sub_in(out, len, acc[len..2 * len], len);\n        return -1;"),
    ("unit-exponent", "bigint.round_bits(out, limbs(k), 0 - 1074 * k, 0);", "bigint.round_bits(out, limbs(k), 0 - 1074, 0);"),
    ("overflow-code", "    if code != 0 {\n        return (0, 1);\n    }\n    if sign < 0 {", "    if sign < 0 {"),
    ("negative-sign-bit", "return (b | (0 - 9223372036854775807 - 1), 0);", "return (b, 0);"),
    ("zero-mantissa-row", "    if m == 0 {\n        return 0;\n    }\n    bigint.clear", "    bigint.clear"),
]

EXPLAINED = {"zero-mantissa-row": "equivalent: the early return for zero only skips adding the number zero (all work limbs are zero)"}


def main():
    t = (ROOT / F).read_text()
    for n, o, w in M:
        if t.count(o) != 1 or o == w:
            print("mutant %s: pattern occurs %d times" % (n, t.count(o)))
            return 1
    if "--check" in sys.argv:
        print("mutants_psum: %d mutants, every pattern applies once" % len(M)); return 0
    bad, killed = [], 0
    for n, o, w in M:
        with tempfile.TemporaryDirectory() as d:
            d = pathlib.Path(d)
            for x in ("lib", "tools"): shutil.copytree(ROOT / x, d / x)
            shutil.copy(ROOT / "cancho.toml", d / "cancho.toml")
            (d / "build").mkdir(); shutil.copytree(ROOT / "build" / "deps", d / "build" / "deps")
            f = d / F; f.write_text(f.read_text().replace(o, w, 1))
            r = subprocess.run([CANCHO, "build", "--bin", "psum-check"], cwd=d, capture_output=True, text=True)
            if r.returncode != 0:
                bad.append(n); print("%-28s INVALID: %s" % (n, r.stderr[-120:].replace("\n", " "))); continue
            try:
                g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_psum.py"), "--driver", str(d / "build" / "psum-check"), "--vectors", "150"], capture_output=True, text=True, timeout=300)
                dead = g.returncode != 0
            except subprocess.TimeoutExpired:
                dead = True
        if dead: killed += 1; print("%-28s killed" % n)
        elif n in EXPLAINED: print("%-28s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-28s SURVIVED" % n)
    print("mutants_psum: %d killed of %d, %d survived and explained, %d unexplained" % (killed, len(M), len(M) - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
