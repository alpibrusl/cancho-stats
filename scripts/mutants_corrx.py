#!/usr/bin/env python3
"""Mutants of `stats corr --method kendall --exact` (docs/rank.md section 15): kendall_p of lib/stest/exact.cho (the counts of the permutations by inversions, the tail, the rounding), the exact parts of put_pair_kendall in lib/stest/answer.cho and of the corr verb in tools/stats/main.cho. Each is applied to a copy, the program rebuilt, scripts/gate_corr_exact.py and the refusals gate rerun; killed = a gate fails or the program traps or hangs. A mutant that does not compile is INVALID. The gates first run on the UNMUTATED build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.
    python3 scripts/mutants_corrx.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
W = "lib/stest/exact.cho"
E = "lib/stest/answer.cho"
R = "lib/csvx/run.cho"
F = "tools/stats/main.cho"
M = [
    ("recurrence-left", W, "if k > 0 {\n                        bigint.add_in(t[co + k * l..", "if k > 1 {\n                        bigint.add_in(t[co + k * l.."),
    ("recurrence-window", W, "if k >= m {\n                        bigint.sub_in(", "if k > m {\n                        bigint.sub_in("),
    ("recurrence-offset", W, "t[po + (k - m) * l..po + (k - m) * l + l], l);", "t[po + (k - m + 1) * l..po + (k - m + 1) * l + l], l);"),
    ("recurrence-start", W, "let t = contents(abw);\n            t[0] = 1;", "let t = contents(abw);\n            t[0] = 2;"),
    ("recurrence-rounds", W, "while m <= n {", "while m < n {"),
    ("discordant", W, "var d = (n0 - s) / 2;", "var d = (n0 - s) / 3;"),
    ("symmetry", W, "if n0 - d < d {", "if n0 - d > d {"),
    ("tail-bound", W, "if k2 <= d {", "if k2 < d {"),
    ("total-bound", W, "while k2 <= n0 {", "while k2 < n0 {"),
    ("doubling", W, "bigint.add_in(sm[2 * l..3 * l], l, sm[0..l], l);\n            bigint.add_in(sm[2 * l..3 * l], l, sm[0..l], l);", "bigint.add_in(sm[2 * l..3 * l], l, sm[0..l], l);"),
    ("cap", W, ">= 0 {\n                pb = 4607182418800017408;", "> 0 {\n                pb = 4607182418800017408;"),
    ("eligible-ties", E, "let exact_data = gx == 0 && gy == 0 && n <= 60;", "let exact_data = gx == 0 || gy == 0 && n <= 60;"),
    ("eligible-limit", E, "let exact_data = gx == 0 && gy == 0 && n <= 60;", "let exact_data = gx == 0 && gy == 0 && n <= 61;"),
    ("auto-uses-exact", E, "if exact_mode >= 1 && exact_data && status == 0 && n >= 3 {", "if exact_mode >= 2 && exact_data && status == 0 && n >= 3 {"),
    ("rows-needed", E, "if exact_mode >= 1 && exact_data && status == 0 && n >= 3 {", "if exact_mode >= 1 && exact_data && status == 0 && n >= 4 {"),
    ("unavailable", E, "let unavailable = exact_mode == 2 && !exact_data;", "let unavailable = exact_mode == 1 && !exact_data;"),
    ("method-name", E, "w = json.put_string(heap, w, \"exact\");", "w = json.put_string(heap, w, \"exactx\");"),
    ("exact-value", E, "p = float_of_bits(xb);", "p = float_of_bits(zb & 9223372036854775807);"),
    ("mode-yes", F, "exact_mode = 2;", "exact_mode = 1;"),
    ("mode-auto", F, "&& cli.has(parsed, table, \"exact\") {\n                                                exact_mode = 1;", "&& !cli.has(parsed, table, \"exact\") {\n                                                exact_mode = 1;"),
    ("refuse-other-methods", F, "\"corr\") && bytes.equal(cli.text(args, parsed, table, \"exact\"), \"yes\") && !bytes.equal(cli.text(args, parsed, table, \"method\"), \"kendall\") {", "\"corr\") && bytes.equal(cli.text(args, parsed, table, \"exact\"), \"yes\") && false {"),
    ("refuse-pairs", F, "if unavailable_pairs > 0 {", "if unavailable_pairs > 1 {"),
    ("count-pairs", F, "unavailable_pairs = unavailable_pairs + 1;", "unavailable_pairs = unavailable_pairs + 0;"),
]

EXPLAINED = {
    "recurrence-start": "equivalent: every count doubles, so does their sum, and 2 lower / total is the same quotient (the counts stay far below the 312 bits of the limbs)",
    "cap": "equivalent: when twice the lower tail equals the total the quotient is exactly 1.0, the value the cap writes",
}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_corr_exact.py"), "--bin", binary, "--files", str(files)], 1500),
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
        print("mutants_corrx: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 100):
        print("mutants_corrx: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_corrx: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
