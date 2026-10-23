#!/usr/bin/env python3
"""Mutants of lib/acc/central.cho (and the three lines of bigint and main that feed it) (docs/describe.md S1-G): each is applied to a copy, the driver rebuilt, scripts/gate_describe.py rerun; killed = the gate fails or the driver traps. A mutant that does
not compile is INVALID. Exit 1 if one survives or is invalid.   python3 scripts/mutants_central.py [--check]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
M = [
    ("k3-n-squared-factor", "lib/acc/central.cho", "limbs_of(n * n, contents(wf));\n                bigint.mul(contents(a1), contents(w3), l3, contents(wf), 3);", "limbs_of(n, contents(wf));\n                bigint.mul(contents(a1), contents(w3), l3, contents(wf), 3);"),
    ("k3-s3-sign", "lib/acc/central.cho", "contents(a1), l3 + 3, s3 < 0);", "contents(a1), l3 + 3, s3 >= 0);"),
    ("k3-three-n-factor", "lib/acc/central.cho", "limbs_of(3 * n, contents(wf));", "limbs_of(2 * n, contents(wf));"),
    ("k3-cross-sign", "lib/acc/central.cho", "contents(a2), l1 + 3 + l2, s1 >= 0);", "contents(a2), l1 + 3 + l2, s1 < 0);"),
    ("k3-cube-factor", "lib/acc/central.cho", "limbs_of(2, contents(wf));\n                bigint.mul(contents(a3), contents(a2), 3 * l1", "limbs_of(3, contents(wf));\n                bigint.mul(contents(a3), contents(a2), 3 * l1"),
    ("k3-cube-sign", "lib/acc/central.cho", "contents(a3), 3 * l1 + 3, s1 < 0);", "contents(a3), 3 * l1 + 3, s1 >= 0);"),
    ("k4-n-cubed", "lib/acc/central.cho", "bigint.mul(contents(a1), contents(wf), 3, contents(wnl), 3);\n                bigint.mul(contents(a2), contents(w4), l4", "bigint.mul(contents(a1), contents(wf), 3, contents(wf), 3);\n                bigint.mul(contents(a2), contents(w4), l4"),
    ("k4-first-term-sign", "lib/acc/central.cho", "contents(a2), l4 + 6, false);", "contents(a2), l4 + 6, true);"),
    ("k4-four-n-squared", "lib/acc/central.cho", "limbs_of(4 * n * n, contents(wf));", "limbs_of(2 * n * n, contents(wf));"),
    ("k4-cross-sign", "lib/acc/central.cho", "s1 < 0 == s3 < 0", "s1 < 0 != s3 < 0"),
    ("k4-six-n", "lib/acc/central.cho", "limbs_of(6 * n, contents(wf));", "limbs_of(5 * n, contents(wf));"),
    ("k4-fourth-power-factor", "lib/acc/central.cho", "limbs_of(3, contents(wf));\n                bigint.mul(contents(a3), contents(a2), 4 * l1", "limbs_of(2, contents(wf));\n                bigint.mul(contents(a3), contents(a2), 4 * l1"),
    ("k4-fourth-power-sign", "lib/acc/central.cho", "contents(a3), 4 * l1 + 3, true);", "contents(a3), 4 * l1 + 3, false);"),
    ("kurtosis-three", "lib/acc/central.cho", "limbs_of(3, contents(wf));\n                bigint.mul(contents(a2), contents(a1), 2 * nl2", "limbs_of(2, contents(wf));\n                bigint.mul(contents(a2), contents(a1), 2 * nl2"),
    ("kurtosis-k4-sign", "lib/acc/central.cho", "mag_len(contents(w41), c), sk4 < 0);", "mag_len(contents(w41), c), sk4 >= 0);"),
    ("kurtosis-subtraction-sign", "lib/acc/central.cho", "contents(a2), 2 * nl2 + 3, true);", "contents(a2), 2 * nl2 + 3, false);"),
    ("l2-subtraction", "lib/acc/central.cho", "contents(a1), 2 * l1, true);", "contents(a1), 2 * l1, false);"),
    ("zero-variance-guard", "lib/acc/central.cho", "if bigint.is_zero(contents(wl), c) {", "if false {"),
    ("skewness-sign", "lib/acc/central.cho", "if sk3 < 0 && c1 == 0 {", "if sk3 > 0 && c1 == 0 {"),
    ("skewness-zero-numerator", "lib/acc/central.cho", "                if sk3 == 0 {\n                    sb = 0;\n                }", ""),
    ("sqrt-too-few-bits", "lib/acc/central.cho", "var k = (132 + bd - bn + 1) / 2;", "var k = (60 + bd - bn + 1) / 2;"),
    ("sqrt-equality", "lib/acc/central.cho", "contents(bw), ns) <= 0 {", "contents(bw), ns) < 0 {"),
    ("sqrt-sticky", "lib/acc/central.cho", "contents(bw), ns) != 0 {\n                            sticky = 1;", "contents(bw), ns) == 0 {\n                            sticky = 1;"),
    ("sqrt-exponent", "lib/acc/central.cho", "round_bits(contents(sw), sl, 0 - k, sticky)", "round_bits(contents(sw), sl, 1 - k, sticky)"),
    ("ratio-too-few-bits", "lib/acc/central.cho", "    var k = 66 + bd - bn;", "    var k = 40 + bd - bn;"),
    ("limbs-of-shift", "lib/acc/central.cho", "    out[1] = v >> 26 & 67108863;", "    out[1] = v >> 25 & 67108863;"),
    ("shl-limb-index", "lib/acc/bigint.cho", "        out[q + i] = v & mask();\n        carry = v >> 26;", "        out[q + i] = v & mask();\n        carry = v >> 25;"),
    ("third-power-accumulator", "tools/stats/main.cho", ", 3, bits_of(x), work);", ", 2, bits_of(x), work);"),
    ("fourth-power-accumulator", "tools/stats/main.cho", ", 4, bits_of(x), work);", ", 3, bits_of(x), work);"),
    ("moments-status-ignored", "tools/stats/main.cho", "if st == 0 && sc == 0 && kc == 0 {", "if sc == 0 && kc == 0 {"),
]

EXPLAINED = {
    "skewness-zero-numerator": "the square root of 0 is 0 bits already; the symmetric columns of gate_describe.py (K3 = 0 exactly, variance not 0) pass with and without the line (equivalent, as far as the gates can show)",
    "sqrt-equality": "dropping an exact hit of a bit leaves s one unit short at 66 or more bits and sets the sticky bit, which changes the rounded double only if the exact root is a tie between two doubles (54 significant bits); a column with such a skewness needs about 2^54 rows, so no gate reaches it: argued, not tested",
    "sqrt-sticky": "the sticky bit decides only when the 66-bit root ends exactly at a half of the last double (a tie); the generated and the fixed columns reached none: argued, not tested",
}


import re


def pattern(o):
    """the text of a mutant as a regex that tolerates the indentation `cancho fmt` gives a continuation line"""
    parts = o.split("\n")
    return r"\n\s*".join(re.escape(p.strip()) if i else re.escape(p) for i, p in enumerate(parts))


def apply(t, o, w):
    return re.sub(pattern(o), lambda m: w, t, count=1)


def occurrences(t, o):
    return len(re.findall(pattern(o), t))


def control():
    """the same two gates on the UNMUTATED build must pass, or every kill below would be a gate that cannot run (a flag it does not know exits non-zero too)"""
    b = str(ROOT / "build" / "stats")
    g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe.py"), "--bin", b, "--files", "80"], capture_output=True, text=True, timeout=600)
    g2 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe_refusals.py"), "--bin", b, "--fuzz", "150"], capture_output=True, text=True, timeout=600)
    return g.returncode == 0 and g2.returncode == 0


def main():
    for n, fl, o, w in M:
        t = (ROOT / fl).read_text()
        if occurrences(t, o) != 1 or o == w:
            print("mutant %s: pattern occurs %d times" % (n, occurrences(t, o)))
            return 1
    if "--check" in sys.argv:
        print("mutants_central: %d mutants, every pattern applies once" % len(M)); return 0
    if not control():
        print("mutants_central: the unmutated build fails the gates: refusing to count kills")
        return 1
    print("control (unmutated) passes both gates")
    bad, killed = [], 0
    for n, fl, o, w in M:
        with tempfile.TemporaryDirectory() as d:
            d = pathlib.Path(d)
            for x in ("lib", "tools", "generated"): shutil.copytree(ROOT / x, d / x)
            shutil.copy(ROOT / "cancho.toml", d / "cancho.toml")
            (d / "build").mkdir(); shutil.copytree(ROOT / "build" / "deps", d / "build" / "deps")
            f = d / fl; f.write_text(apply(f.read_text(), o, w))
            r = subprocess.run([CANCHO, "build", "--bin", "stats"], cwd=d, capture_output=True, text=True)
            if r.returncode != 0:
                bad.append(n); print("%-28s INVALID: %s" % (n, r.stderr[-120:].replace("\n", " "))); continue
            try:
                g = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe.py"), "--bin", str(d / "build" / "stats"), "--files", "80"], capture_output=True, text=True, timeout=300)
                dead = g.returncode != 0
                if not dead:
                    g2 = subprocess.run([sys.executable, str(ROOT / "scripts" / "gate_describe_refusals.py"), "--bin", str(d / "build" / "stats"), "--fuzz", "150"], capture_output=True, text=True, timeout=600)
                    dead = g2.returncode != 0
            except subprocess.TimeoutExpired:
                dead = True
        if dead: killed += 1; print("%-28s killed" % n)
        elif n in EXPLAINED: print("%-28s SURVIVED (explained: %s)" % (n, EXPLAINED[n]))
        else: bad.append(n); print("%-28s SURVIVED" % n)
    print("mutants_central: %d killed of %d, %d survived and explained, %d unexplained" % (killed, len(M), len(M) - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
