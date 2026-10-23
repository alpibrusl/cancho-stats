#!/usr/bin/env python3
"""Mutants of lib/acc/xsum.cho (docs/corr.md C-I): each is applied to a copy, the driver rebuilt, scripts/gate_xsum.py rerun; killed = the gate fails or the driver traps. A mutant that does
not compile is INVALID. Exit 1 if one survives or is invalid.   python3 scripts/mutants_xsum.py [--check]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
F = "lib/acc/xsum.cho"
M = [
    ("hidden-bit-dropped-x", "    if ef != 0 {\n        m = m | 4503599627370496;\n        f = ef - 1;\n    }\n    return (m, f);", "    if ef != 0 {\n        f = ef - 1;\n    }\n    return (m, f);"),
    ("exponent-offset", "        f = ef - 1;\n    }\n    return (m, f);", "        f = ef;\n    }\n    return (m, f);"),
    ("shift-ignores-y", "return psum.add_shifted(acc, base, len, work, 6, 6, fx + fy);", "return psum.add_shifted(acc, base, len, work, 6, 6, fx + fx);"),
    ("shift-ignores-x", "return psum.add_shifted(acc, base, len, work, 6, 6, fx + fy);", "return psum.add_shifted(acc, base, len, work, 6, 6, fy + fy);"),
    ("sign-is-x-only", "if (bx < 0) != (by < 0) {", "if bx < 0 {"),
    ("sign-is-and", "if (bx < 0) != (by < 0) {", "if bx < 0 && by < 0 {"),
    ("sign-is-or", "if (bx < 0) != (by < 0) {", "if bx < 0 || by < 0 {"),
    ("sign-inverted", "if (bx < 0) != (by < 0) {", "if (bx < 0) == (by < 0) {"),
    ("mantissa-limb-1", "    work[1] = mx >> 26 & mask();", "    work[1] = mx >> 27 & mask();"),
    ("mantissa-limb-2", "    work[2] = mx >> 52;", "    work[2] = mx >> 51;"),
    ("y-mantissa-limb-1", "    work[4] = my >> 26 & mask();", "    work[4] = my >> 27 & mask();"),
    ("y-mantissa-limb-2", "    work[5] = my >> 52;", "    work[5] = my >> 51;"),
    ("y-mantissa-limb-0", "    work[3] = my & mask();", "    work[3] = mx & mask();"),
    ("product-operands", "bigint.mul(work[6..12], work[0..3], 3, work[3..6], 3);", "bigint.mul(work[6..12], work[0..3], 3, work[0..3], 3);"),
    ("product-short", "bigint.mul(work[6..12], work[0..3], 3, work[3..6], 3);", "bigint.mul(work[6..12], work[0..3], 3, work[3..6], 2);"),
    ("product-offset", "return psum.add_shifted(acc, base, len, work, 6, 6, fx + fy);", "return psum.add_shifted(acc, base, len, work, 6, 5, fx + fy);"),
    ("zero-x-row", "    if mx == 0 || my == 0 {\n        return 0;\n    }", "    if my == 0 {\n        return 0;\n    }"),
    ("zero-y-row", "    if mx == 0 || my == 0 {\n        return 0;\n    }", "    if mx == 0 {\n        return 0;\n    }"),
    ("size-of-k1", "return psum.size(2);", "return psum.size(1);"),
    ("rounded-power", "return psum.rounded(acc, 2, out);", "return psum.rounded(acc, 1, out);"),
]

EXPLAINED = {
    "product-offset": "equivalent: the product of two 53-bit mantissas is below 2^106, so the sixth 26-bit limb (bits 130 and up) is always zero and adding five limbs adds the same number",
    "zero-x-row": "equivalent: with mx = 0 the product is zero and adding the zero limbs changes nothing; the early return only saves the work",
    "zero-y-row": "equivalent: as zero-x-row, for my = 0",
}


def control():
    """the gate on the UNMUTATED build must pass, or every kill below would be a gate that cannot run"""
    g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_xsum.py"), "--vectors", "150"], capture_output=True, text=True, timeout=600)
    return g.returncode == 0


def main():
    t = (ROOT / F).read_text()
    for n, o, w in M:
        if t.count(o) != 1 or o == w:
            print("mutant %s: pattern occurs %d times" % (n, t.count(o)))
            return 1
    if "--check" in sys.argv:
        print("mutants_xsum: %d mutants, every pattern applies once" % len(M)); return 0
    if not control():
        print("mutants_xsum: the unmutated build fails the gate: refusing to count kills")
        return 1
    print("control (unmutated) passes the gate")
    bad, killed = [], 0
    for n, o, w in M:
        with tempfile.TemporaryDirectory() as d:
            d = pathlib.Path(d)
            for x in ("lib", "tools"): shutil.copytree(ROOT / x, d / x)
            shutil.copy(ROOT / "cancho.toml", d / "cancho.toml")
            (d / "build").mkdir(); shutil.copytree(ROOT / "build" / "deps", d / "build" / "deps")
            f = d / F; f.write_text(f.read_text().replace(o, w, 1))
            r = subprocess.run([CANCHO, "build", "--bin", "xsum-check"], cwd=d, capture_output=True, text=True)
            if r.returncode != 0:
                bad.append(n); print("%-28s INVALID: %s" % (n, r.stderr[-120:].replace("\n", " "))); continue
            try:
                g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_xsum.py"), "--driver", str(d / "build" / "xsum-check"), "--psum", str(ROOT / "build" / "psum-check"), "--vectors", "150"], capture_output=True, text=True, timeout=300)
                dead = g.returncode != 0
            except subprocess.TimeoutExpired:
                dead = True
        if dead: killed += 1; print("%-28s killed" % n)
        elif n in EXPLAINED: print("%-28s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-28s SURVIVED" % n)
    print("mutants_xsum: %d killed of %d, %d survived and explained, %d unexplained" % (killed, len(M), len(M) - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
