#!/usr/bin/env python3
"""Mutants of `stats test kruskal` (docs/rank.md R-H): lib/stest/kruskal.cho (the ranks, the ties, the exact integers), lib/stest/answer6.cho (the entries), the kruskal parts of lib/csvx/run.cho and of tools/stats/main.cho.
Each is applied to a copy, the program rebuilt, scripts/gate_kruskal.py and the refusals gate rerun; killed = a gate fails or the program traps or hangs. A mutant that does not compile is INVALID. The gates first run on the
UNMUTATED build and the script refuses to count kills if they fail. Exit 1 if one survives unexplained or is invalid.   python3 scripts/mutants_kruskal.py [--check] [--only NAME,NAME]"""
import os, pathlib, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CANCHO = os.environ.get("CANCHO", str(ROOT.parent / "cancho-compiler-a4572ea" / "target" / "release" / "cancho"))
K = "lib/stest/kruskal.cho"
E = "lib/stest/answer6.cho"
R = "lib/csvx/run.cho"
F = "tools/stats/main.cho"
M = [
    # the ranks and the ties
    ("rank-below", K, "let lo = ranks.lower(contents(sw), m, keys[i]);", "let lo = ranks.upper(contents(sw), m, keys[i]);"),
    ("rank-at-most", K, "let hi = ranks.upper(contents(sw), m, keys[i]);", "let hi = ranks.lower(contents(sw), m, keys[i]);"),
    ("rank-doubled", K, "ds[gid[i]] = ds[gid[i]] + lo + hi + 1;", "ds[gid[i]] = ds[gid[i]] + lo + hi;"),
    ("rank-group", K, "ds[gid[i]] = ds[gid[i]] + lo + hi + 1;", "ds[gid[i]] = ds[gid[0]] + lo + hi + 1;"),
    ("count-group", K, "ns[gid[i]] = ns[gid[i]] + 1;", "ns[gid[i]] = ns[gid[i]] + 2;"),
    ("tie-start", K, "var end = at + 1;", "var end = at + 2;"),
    ("tie-single-excluded", K, "if tt > 1 {\n                        ties = ties + 1;", "if tt > 0 {\n                        ties = ties + 1;"),
    ("tie-size-sum", K, "sum_t = sum_t + tt;", "sum_t = sum_t + 1;"),
    ("tie-square", K, "welch.smul(f, 1, 0, 3, 0, 3);\n                        welch.smul(f, 2, 1, 6, 0, 3);", "welch.smul(f, 1, 0, 3, 0, 3);\n                        welch.smul(f, 2, 1, 6, 1, 3);"),
    ("tie-cube-into", K, "bigint.add_in(t[6 * c..7 * c], c, f[2 * welch.slot()..2 * welch.slot() + 9], 9);", "bigint.add_in(t[6 * c..7 * c], c, f[welch.slot()..welch.slot() + 6], 6);"),
    # the fraction and the exact integers
    ("fold-group-count", K, "central.limbs_of(ns[i], f[0..3]);", "central.limbs_of(ns[i] + 1, f[0..3]);"),
    ("fold-rank-sum", K, "central.limbs_of(ds[i], f[welch.slot()..welch.slot() + 3]);", "central.limbs_of(ds[i] / 2, f[welch.slot()..welch.slot() + 3]);"),
    ("fold-square", K, "welch.smul(f, 2, 1, 3, 1, 3);", "welch.smul(f, 2, 1, 3, 0, 3);"),
    ("fold-numerator", K, "welch.mulf(t, 2, 0, welch.len_of(t, 0), f, 0, 3);", "welch.mulf(t, 2, 0, welch.len_of(t, 0), f, 1, 3);"),
    ("fold-denominator", K, "welch.mulf(t, 3, 1, welch.len_of(t, 1), f, 2, 6);", "welch.mulf(t, 3, 0, welch.len_of(t, 0), f, 2, 6);"),
    ("fold-den-update", K, "welch.mulf(t, 2, 1, welch.len_of(t, 1), f, 0, 3);\n                            bigint.copy(t[c..2 * c]", "welch.mulf(t, 2, 1, welch.len_of(t, 1), f, 1, 3);\n                            bigint.copy(t[c..2 * c]"),
    ("groups-needed", K, "if k < 2 {", "if k < 1 {"),
    ("n-plus-one", K, "central.limbs_of(nn + 1, f[4 * welch.slot()..4 * welch.slot() + 3]);", "central.limbs_of(nn, f[4 * welch.slot()..4 * welch.slot() + 3]);"),
    ("n-plus-one-squared", K, "welch.smul(f, 5, 4, 3, 4, 3);", "welch.smul(f, 5, 4, 3, 3, 3);"),
    ("n-cube", K, "welch.smul(f, 8, 7, 6, 3, 3);", "welch.smul(f, 8, 7, 6, 4, 3);"),
    ("x-bound", K, "welch.mulf(t, 5, 1, welch.len_of(t, 1), f, 6, 9);", "welch.mulf(t, 5, 1, welch.len_of(t, 1), f, 5, 6);"),
    ("e-cube-plus-t", K, "bigint.add_in(t[8 * c..9 * c], c, f[11 * welch.slot()..11 * welch.slot() + 3], 3);", "bigint.add_in(t[8 * c..9 * c], c, f[3 * welch.slot()..3 * welch.slot() + 3], 3);"),
    ("e-n-minus", K, "bigint.add_in(t[9 * c..10 * c], c, f[3 * welch.slot()..3 * welch.slot() + 3], 3);", "bigint.add_in(t[9 * c..10 * c], c, f[11 * welch.slot()..11 * welch.slot() + 3], 3);"),
    ("e-zero", K, "if bigint.is_zero(t[8 * c..9 * c], c) {\n                            status = 2;", "if bigint.is_zero(t[7 * c..8 * c], c) {\n                            status = 2;"),
    ("h-factor", K, "central.limbs_of(3 * (nn - 1), f[9 * welch.slot()..9 * welch.slot() + 3]);", "central.limbs_of(3 * nn, f[9 * welch.slot()..9 * welch.slot() + 3]);"),
    ("h-denominator", K, "welch.mulr(t, 10, 1, welch.len_of(t, 1), 8, welch.len_of(t, 8));", "welch.mulr(t, 10, 1, welch.len_of(t, 1), 6, welch.len_of(t, 6));"),
    ("eps-factor", K, "central.limbs_of(3, f[10 * welch.slot()..10 * welch.slot() + 3]);", "central.limbs_of(6, f[10 * welch.slot()..10 * welch.slot() + 3]);"),
    ("h-uses-eps-factor", K, "welch.mulf(t, 11, 7, welch.len_of(t, 7), f, 9, 3);", "welch.mulf(t, 11, 7, welch.len_of(t, 7), f, 10, 3);"),
    ("sum-t-limbs", K, "central.limbs_of(sum_t, f[11 * welch.slot()..11 * welch.slot() + 3]);", "central.limbs_of(sum_t + 1, f[11 * welch.slot()..11 * welch.slot() + 3]);"),
    # the entries
    ("entry-status-one", E, "    if r.status == 1 {\n        reason = \"too-few-rows\";", "    if r.status == 2 {\n        reason = \"too-few-rows\";"),
    ("entry-status-two", E, "} else if r.status == 2 {\n        reason = \"variance-zero\";", "} else if r.status == 1 {\n        reason = \"variance-zero\";"),
    ("entry-df", E, "dist.chi2_sf(float_of_bits(r.hb), float_of(r.k - 1))", "dist.chi2_sf(float_of_bits(r.hb), float_of(r.k))"),
    ("entry-df-printed", E, "w = json.put_int(heap, w, r.k - 1);", "w = json.put_int(heap, w, r.k);"),
    ("entry-df-null", E, "    if r.status != 1 {\n        w = json.put_int", "    if r.status != 2 {\n        w = json.put_int"),
    ("entry-p-method", E, "w = json.put_string(heap, w, \"chi2\");", "w = json.put_string(heap, w, \"f\");"),
    ("entry-effect-name", E, "w = json.put_string(heap, w, \"epsilon_squared\");", "w = json.put_string(heap, w, \"eps\");"),
    ("entry-effect-value", E, "w = json.put_float(heap, w, float_of_bits(r.eb));", "w = json.put_float(heap, w, float_of_bits(r.hb));"),
    ("entry-mean-rank", E, "float_of(ds[i]) / (2.0 * float_of(ns[i]))", "float_of(ds[i]) / float_of(ns[i])"),
    ("entry-mean-rank-empty", E, "if ns[i] > 0 {\n            w = json.put_float", "if ns[i] > 1 {\n            w = json.put_float"),
    ("entry-ties", E, "w = json.put_int(heap, w, r.ties);", "w = json.put_int(heap, w, r.k);"),
    ("entry-method", E, "w = json.put_string(heap, w, \"kruskal\");", "w = json.put_string(heap, w, \"anova\");"),
    ("entry-beyond", E, "} else if r.hc != 0 {", "} else if r.hc != 1 {"),
    # the driver
    ("run-kind-name", R, "w = json.put_string(heap, w, \"kruskal\");", "w = json.put_string(heap, w, \"anova\");"),
    ("run-collect-group", R, "gi = vec.push(heap, gi, i);", "gi = vec.push(heap, gi, 0);"),
    ("run-collect-column", R, "gk = vec.push(heap, gk, vec.get(store, r * nc + j));", "gk = vec.push(heap, gk, vec.get(store, r * nc));"),
    ("run-present", R, "if key != rank.min_int() {\n                                    contents(pkw)[m] = key;", "if key == rank.min_int() {\n                                    contents(pkw)[m] = key;"),
    ("run-stride", R, "                                let key = vec.get(gk, r * nsel + si);", "                                let key = vec.get(gk, r * nsel);"),
    ("run-cell-limit", R, "if held > cli.nat(args, parsed, table, \"max-cells\") {", "if held >= cli.nat(args, parsed, table, \"max-cells\") {"),
    ("run-numeric-total", R, "contents(tnw)[jj] = contents(tnw)[jj] + contents(cv.numeric)[jj];", "contents(tnw)[jj] = contents(tnw)[jj] + 0;"),
    ("run-bad-total", R, "contents(tbw)[jj] = contents(tbw)[jj] + contents(cv.bad)[jj];", "contents(tbw)[jj] = contents(tbw)[jj] + 0;"),
    ("main-group-required", F, "\"kruskal\") && !cli.has(parsed, table, \"group\") {", "\"kruskal\") && false {"),
    ("main-dispatch", F, "if bytes.equal(cli.operand(args, parsed, 1), \"kruskal\") {\n                            kind_of_groups = 1;", "if bytes.equal(cli.operand(args, parsed, 1), \"kruskalx\") {\n                            kind_of_groups = 1;"),
]

EXPLAINED = {}


def gates(binary, files):
    """the gates in the order that fails fastest; True if all pass"""
    for argv, tmo in (([sys.executable, str(ROOT / "scripts" / "gate_kruskal.py"), "--bin", binary, "--files", str(files)], 1200),
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
        print("mutants_kruskal: %d mutants, every pattern applies once" % len(M)); return 0
    if not gates(str(ROOT / "build" / "stats"), 100):
        print("mutants_kruskal: the unmutated build fails the gates: refusing to count kills")
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
    print("mutants_kruskal: %d killed of %d, %d survived and explained, %d unexplained" % (killed, total, total - killed - len(bad), len(bad)))
    if bad: print("SURVIVORS: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
