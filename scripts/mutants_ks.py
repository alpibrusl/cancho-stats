#!/usr/bin/env python3
"""Mutants of `stats test ks` (docs/rank.md R-H): lib/stest/ks.cho (the walk, D, lambda, the exact p), lib/stest/exact.cho (the path count), lib/stest/answer9.cho (the entries), the ks parts of
lib/csvx/run.cho and of tools/stats/main.cho. Each is applied to a copy, the program rebuilt, scripts/gate_ks.py and the refusals gate rerun; killed = a gate fails or the program traps or hangs. A mutant that does
not compile is INVALID. The gates first run on the UNMUTATED build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.
    python3 scripts/mutants_ks.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
W = "lib/stest/ks.cho"
X = "lib/stest/exact.cho"
E = "lib/stest/answer9.cho"
R = "lib/csvx/run.cho"
F = "tools/stats/main.cho"
M = [
    # the walk and D
    ("d-sign", W, "var dev = i * nb - j * na;", "var dev = i * na - j * nb;"),
    ("d-abs", W, "if dev < 0 {\n                                dev = 0 - dev;", "if dev > 0 {\n                                dev = 0 - dev;"),
    ("d-max", W, "if dev > d {", "if dev < d {"),
    ("tie-step", W, "while end < m && keys[o[end]] == keys[o[at]] {", "while end < m && keys[o[end]] >= keys[o[at]] {"),
    ("ties-count", W, "if end - at > 1 {", "if end - at > 2 {"),
    ("group-count", W, "if gid[r] == 0 {\n                            na = na + 1;", "if gid[r] != 0 {\n                            na = na + 1;"),
    ("group-walk", W, "if gid[o[q]] == 0 {\n                                    i = i + 1;", "if gid[o[q]] != 0 {\n                                    i = i + 1;"),
    ("groups-needed", W, "if na == 0 || nb == 0 {", "if na == 0 && nb == 0 {"),
    ("dstat", W, "dstat = float_of(d) / float_of(na * nb);", "dstat = float_of(d) / float_of(na + nb);"),
    # lambda
    ("lambda-n", W, "central.limbs_of(m, f[2 * welch.slot()..2 * welch.slot() + 3]);", "central.limbs_of(m - 1, f[2 * welch.slot()..2 * welch.slot() + 3]);"),
    ("lambda-product", W, "central.limbs_of(na * nb, f[3 * welch.slot()..3 * welch.slot() + 3]);", "central.limbs_of(na + nb, f[3 * welch.slot()..3 * welch.slot() + 3]);"),
    ("lambda-square", W, "welch.smul(f, 1, 0, 3, 0, 3);", "welch.smul(f, 1, 0, 3, 3, 3);"),
    # the exact p
    ("mode-auto", W, "use_exact = mannwhitney.possible(keys, m);", "use_exact = false;"),
    ("mode-yes", W, "} else if mode == 1 {\n                            use_exact = true;", "} else if mode == 2 {\n                            use_exact = true;"),
    ("d-zero", W, "if d == 0 {\n                                pb = bits_of(1.0);", "if d == 1 {\n                                pb = bits_of(1.0);"),
    ("p-complement", W, "central.limbs_of(total - inside, f[5 * welch.slot()..5 * welch.slot() + 3]);", "central.limbs_of(inside, f[5 * welch.slot()..5 * welch.slot() + 3]);"),
    ("paths-band-start", X, "if dev < d {\n                            a[0] = 1;", "if dev <= d {\n                            a[0] = 1;"),
    ("paths-band", X, "if dev < d {\n                            a[i * width + j] = left + up;", "if dev <= d {\n                            a[i * width + j] = left + up;"),
    ("paths-total", X, "b[i * width + j] = leftb + upb;", "b[i * width + j] = leftb;"),
    ("paths-total-step", X, "leftb = b[(i - 1) * width + j];", "leftb = b[(i - 1) * width + j] * 2;"),
    ("paths-start", X, "b[0] = 1;", "b[0] = 2;"),
    ("paths-end", X, "inside = a[na * width + nb];", "inside = a[na * width];"),
    ("paths-up", X, "up = a[i * width + j - 1];", "up = a[i * width + j - 1] + a[i * width + j - 1];"),
    # the entries
    ("entry-status", E, "    if r.status == 1 {\n        reason = \"too-few-rows\";", "    if r.status == 2 {\n        reason = \"too-few-rows\";"),
    ("entry-reason", E, "reason = \"too-few-rows\";", "reason = \"p-not-computed\";"),
    ("entry-exact", E, "} else if r.exact == 1 {", "} else if r.exact == 0 {"),
    ("entry-method", E, "method = \"exact\";", "method = \"asymptotic\";"),
    ("entry-pc", E, "if r.pc == 0 {", "if r.pc == 1 {"),
    ("entry-asymptotic", E, "} else if r.lc == 0 {", "} else if r.lc == 1 {"),
    ("entry-code", E, "if code == 0 {\n            have_p = true;\n            p = q;", "if code == 1 {\n            have_p = true;\n            p = q;"),
    ("entry-statistic", E, "if r.status != 1 {\n        w = json.put_float(heap, w, r.dstat);", "if r.status != 2 {\n        w = json.put_float(heap, w, r.dstat);"),
    ("entry-distance", E, "w = json.put_int(heap, w, r.d);", "w = json.put_int(heap, w, r.na);"),
    ("entry-n", E, "w = json.put_int(heap, w, r.na);\n    w = json.put_key(heap, w, \"b\");", "w = json.put_int(heap, w, r.nb);\n    w = json.put_key(heap, w, \"b\");"),
    ("entry-ties", E, "w = json.put_int(heap, w, r.ties);", "w = json.put_int(heap, w, r.d);"),
    ("entry-method-name", E, "w = json.put_string(heap, w, \"ks\");", "w = json.put_string(heap, w, \"ksx\");"),
    ("entry-lambda", E, "w = json.put_float(heap, w, float_of_bits(r.lb));\n    } else {\n        w = json.put_null(heap, w);\n    }\n    w = json.put_key(heap, w, \"p_value\");", "w = json.put_float(heap, w, r.dstat);\n    } else {\n        w = json.put_null(heap, w);\n    }\n    w = json.put_key(heap, w, \"p_value\");"),
    # the driver
    ("run-levels", R, "if clean && kind >= 2 && k != 2 {", "if clean && kind == 2 && k != 2 {"),
    ("run-precheck", R, "if clean && kind >= 2 && bytes.equal(cli.text(args, parsed, table, \"exact\"), \"yes\") {", "if clean && kind == 2 && bytes.equal(cli.text(args, parsed, table, \"exact\"), \"yes\") {"),
    ("run-kind-name", R, "w = json.put_string(heap, w, \"ks\");", "w = json.put_string(heap, w, \"mannwhitney\");"),
    ("run-dispatch", R, "if kind == 3 {\n                                borrow sel as &sr8 in {", "if kind == 4 {\n                                borrow sel as &sr8 in {"),
    ("run-mode-yes", R, "if bytes.equal(exact_text, \"yes\") {\n        mode = 1;\n    } else if bytes.equal(exact_text, \"no\") {\n        mode = 2;\n    }\n    var si = 0;", "if bytes.equal(exact_text, \"no\") {\n        mode = 1;\n    } else if bytes.equal(exact_text, \"no\") {\n        mode = 2;\n    }\n    var si = 0;"),
    ("run-mode-no", R, "} else if bytes.equal(exact_text, \"no\") {\n        mode = 2;\n    }\n    var si = 0;", "} else if bytes.equal(exact_text, \"yes\") {\n        mode = 2;\n    }\n    var si = 0;"),
    ("run-mode", R, "let kr = ks.run(heap, contents(pkw), contents(pgw), m, mode);", "let kr = ks.run(heap, contents(pkw), contents(pgw), m, 0);"),
    ("run-count", R, "let kr = ks.run(heap, contents(pkw), contents(pgw), m, mode);", "let kr = ks.run(heap, contents(pkw), contents(pgw), m - 1, mode);"),
    ("run-stride", R, "let key = vec.get(gk, r * nsel + si);\n                        if key != rank.min_int() {\n                            contents(pkw)[m] = key;\n                            contents(pgw)[m] = vec.get(gi, r);\n                            m = m + 1;\n                        }\n                        r = r + 1;\n                    }\n                    let kr = ks.run", "let key = vec.get(gk, r * nsel);\n                        if key != rank.min_int() {\n                            contents(pkw)[m] = key;\n                            contents(pgw)[m] = vec.get(gi, r);\n                            m = m + 1;\n                        }\n                        r = r + 1;\n                    }\n                    let kr = ks.run"),
    ("main-group-required", F, "bytes.equal(cli.operand(args, parsed, 1), \"ks\") && !cli.has(parsed, table, \"group\") {", "bytes.equal(cli.operand(args, parsed, 1), \"ks\") && false {"),
    ("main-dispatch", F, "} else if bytes.equal(cli.operand(args, parsed, 1), \"ks\") {\n                            kind_of_groups = 3;", "} else if bytes.equal(cli.operand(args, parsed, 1), \"ksx\") {\n                            kind_of_groups = 3;"),
]

EXPLAINED = {
}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_ks.py"), "--bin", binary, "--files", str(files)], 1500),
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
        print("mutants_ks: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 100):
        print("mutants_ks: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_ks: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
