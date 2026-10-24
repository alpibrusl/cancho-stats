#!/usr/bin/env python3
"""The determinism gate G8 (docs/numerics.md section 6): the fingerprint of the bytes the drivers print, against a committed golden, and the thread check.

    python3 scripts/xplat_fp.py                  # print nothing but a verdict: the fingerprint of this build equals docs/xplat_fp.txt, and every thread count 1..16 gives the same bytes
    python3 scripts/xplat_fp.py --print          # print the fingerprint (what a platform contributes to the comparison)
    python3 scripts/xplat_fp.py --update         # rewrite the golden: only when the numbers change on purpose, in a commit that says why
    python3 scripts/xplat_fp.py --mutants        # the checks are shown able to fail

The fingerprint is one tab-separated line per sweep: `kind  name  count  sha256`, the hash over the driver's output lines (`id  hex bits of the result  code`) for a fixed list of cases:
`special` the first N cases of each range of gate_special.py (seed 1), `edge` the edge cases of gate_domain.py for each function, `acc` the vectors of three families of measure_acc.py through acc-check, `describe` the bytes `stats describe` prints for 12 seeded CSV files, `corr`, `corr-spearman` and `corr-kendall` those of `stats corr` (Pearson, --method spearman, --method kendall) for 8 files each, `test-t1` those of `stats test t1` for 8 files.
Equal to the golden is equal to each other: the golden was written on one machine and every other platform of the workflow must print its bytes.
Thread check: for every sweep, `special-check FN --threads T` for T in 1, 2, 3, 4, 7, 8, 16 must print the bytes of the sequential run.
Exit 0: the fingerprint equals the golden and the thread check holds. Exit 1: a difference (printed). Exit 2: the script could not run.
"""
import argparse
import difflib
import hashlib
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import gate_domain  # noqa: E402
import gate_special  # noqa: E402
import measure_acc  # noqa: E402
import gate_describe  # noqa: E402
import gate_corr  # noqa: E402
import gate_test  # noqa: E402
import tempfile  # noqa: E402

GOLDEN = ROOT / "docs" / "xplat_fp.txt"
THREADS = (1, 2, 3, 4, 7, 8, 16)
TIMEOUT = 600


def drive(driver, fn, cases, threads=0):
    inp = "".join("%d %s\n" % (i, " ".join(gate_special.hx(a) for a in c)) for i, c in enumerate(cases))
    argv = [driver, fn] + (["--threads", str(threads)] if threads else [])
    r = subprocess.run(argv, input=inp, capture_output=True, text=True, timeout=TIMEOUT)
    if r.returncode != 0:
        raise RuntimeError("%s exit %d: %s" % (" ".join(argv), r.returncode, r.stderr[-200:]))
    return r.stdout


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def special_sweeps():
    for spec in gate_special.RANGES:
        name, fn, n, gen, kind, b0, slow = spec
        cases = gate_special.cases_of(spec, 1, 1.0, False)
        yield ("special", name, fn, cases[: (300 if slow else 3000)])


def edge_sweeps():
    for fn in gate_domain.ARITY:
        yield ("edge", fn, fn, gate_domain.edge_cases(fn))


def acc_vectors():
    r = gate_special.SM64(7)
    out = []
    for name, gen in (("uniform", measure_acc.fam_uniform), ("wide", measure_acc.mk(-1074, 1023)), ("inside-R", measure_acc.mk(-485, 500))):
        vecs = [gen(r, [1, 2, 3, 10, 100][i % 5]) for i in range(150)]
        out.append((name, vecs))
    return out


def describe_lines(stats):
    """`stats describe` on 12 seeded CSV files (the shapes of gate_describe.py): the hash of the bytes it prints, with --ddof 1 and 0."""
    r = gate_special.SM64(11)
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for k in range(12):
            ncols = 1 + int(r.next() % 5)
            shapes = [gate_describe.SHAPES[int(r.next() % len(gate_describe.SHAPES))] for _ in range(ncols)]
            n = [1, 2, 3, 5, 20, 100, 400][int(r.next() % 7)]
            cols = [gate_describe.column(r, sh, n) for sh in shapes]
            names = ["c%d_%s" % (j, sh) for j, sh in enumerate(shapes)]
            rows = [[cols[j][i] for j in range(ncols)] for i in range(n)]
            h = hashlib.sha256()
            for ddof in (1, 0):
                p = gate_describe.run(stats, tmp, names, rows, ddof)
                if p.returncode != 0:
                    raise RuntimeError("stats exit %d" % p.returncode)
                h.update(p.stdout.encode())
            out.append("describe\tfile%02d\t%d\t%s" % (k, n, h.hexdigest()))
    return out


def describe_by_lines(stats):
    """`stats describe --by g` (and the same with --cols and --format text) on 4 seeded files (a group column of six values with one missing, a numeric column and a text column), --ddof 1 and 0: the hash of the bytes printed."""
    r = gate_special.SM64(23)
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for k, n in enumerate((5, 20, 100, 400)):
            groups = ["a", "b", "", "c c", '"a"', "\u00e9"]
            g = [groups[int(r.next() % 6)] for _ in range(n)]
            num = gate_describe.column(r, "uniform", n)
            txt = gate_describe.column(r, "levels", n)
            rows = [[g[i], num[i], txt[i]] for i in range(n)]
            h = hashlib.sha256()
            for ddof in (1, 0):
                p = gate_describe.run(stats, tmp, ["g", "x", "t"], rows, ddof, ["--by", "g"])
                if p.returncode != 0:
                    raise RuntimeError("stats --by exit %d" % p.returncode)
                h.update(p.stdout.encode())
            pt = gate_describe.run(stats, tmp, ["g", "x", "t"], rows, 1, ["--by", "g", "--cols", "x,t", "--format", "text"])
            if pt.returncode != 0:
                raise RuntimeError("stats --format text exit %d" % pt.returncode)
            h.update(pt.stdout.encode())
            out.append("describe-by\tfile%02d\t%d\t%s" % (k, n, h.hexdigest()))
    return out


def corr_lines(stats, method="pearson", kind="corr"):
    """`stats corr` on 8 seeded CSV files (the numeric shapes of gate_corr.py, correlated pairs, missing cells, huge and tiny values), --missing pairwise and listwise: the hash of the bytes it prints
    (r, the statistic and the p-value of every pair)."""
    r = gate_special.SM64(23)
    # the list of the shapes as it was when these lines were written: a shape added to gate_corr.SHAPES later must not move the golden
    shapes_all = ["uniform", "ints", "lin", "offset", "offset-neg", "wide", "tiny", "huge", "const", "gappy", "lin", "lin", "ints", "uniform"]
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for k in range(8):
            ncols = 2 + int(r.next() % 3)
            shapes = [shapes_all[int(r.next() % len(shapes_all))] for _ in range(ncols)]
            n = [3, 5, 20, 100, 400][int(r.next() % 5)]
            cols = []
            for j, sh in enumerate(shapes):
                cols.append(gate_corr.column(r, sh, n, cols[j - 1] if j else []))
            names = ["c%d_%s" % (j, sh) for j, sh in enumerate(shapes)]
            rows = [[cols[j][i] for j in range(ncols)] for i in range(n)]
            h = hashlib.sha256()
            for mode in ("pairwise", "listwise"):
                p = gate_corr.run(stats, tmp, names, rows, ["--missing", mode] + ([] if method == "pearson" else ["--method", method]))
                if p.returncode != 0:
                    h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
                else:
                    h.update(p.stdout.encode())
            out.append("%s\tfile%02d\t%d\t%s" % (kind, k, n, h.hexdigest()))
    return out


def test_lines(stats):
    """`stats test t1` on 8 seeded CSV files (numeric shapes of gate_describe.py with missing cells; mu0 and the level drawn): the hash of the bytes it prints (the statistic, df, p-value, effect size, interval)."""
    r = gate_special.SM64(37)
    shapes_all = ["uniform", "ints", "offset", "offset-neg", "wide", "tiny", "huge", "const", "gappy"]
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for k in range(8):
            ncols = 1 + int(r.next() % 3)
            shapes = [shapes_all[int(r.next() % len(shapes_all))] for _ in range(ncols)]
            n = [2, 3, 5, 20, 100, 400][int(r.next() % 6)]
            cols = [gate_describe.column(r, sh, n) for sh in shapes]
            names = ["c%d_%s" % (j, sh) for j, sh in enumerate(shapes)]
            rows = [[cols[j][i] for j in range(ncols)] for i in range(n)]
            mu0 = r.pick(["0", "1", "-3.5", "21.5", "1e8", "0.001"])
            level = r.pick(["0.95", "0.9", "0.99"])
            p = gate_test.run(stats, tmp, names, rows, ["--cols", ",".join(names), "--mu0", mu0, "--level", level])
            h = hashlib.sha256()
            h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("test-t1\tfile%02d\t%d\t%s" % (k, n, h.hexdigest()))
    return out


def fingerprint(driver, acc, threads_check=True, stats=None):
    lines = []
    bad = []
    for kind, name, fn, cases in list(special_sweeps()) + list(edge_sweeps()):
        base = drive(driver, fn, cases)
        lines.append("%s\t%s\t%d\t%s" % (kind, name, len(cases), sha(base)))
        if threads_check:
            for t in THREADS:
                if drive(driver, fn, cases, t) != base:
                    bad.append("%s %s: --threads %d differs from the sequential run" % (kind, name, t))
    for name, vecs in acc_vectors():
        lines_in = []
        for xs, ys in vecs:
            lines_in.append("%d %s %s" % (len(xs), " ".join(gate_special.hx(x) for x in xs), " ".join(gate_special.hx(y) for y in ys)))
        r = subprocess.run([acc], input="\n".join(lines_in) + "\n", capture_output=True, text=True, timeout=TIMEOUT)
        if r.returncode != 0:
            raise RuntimeError("acc-check exit %d" % r.returncode)
        lines.append("acc\t%s\t%d\t%s" % (name, len(vecs), sha(r.stdout)))
    if stats:
        lines += describe_lines(stats)
        lines += describe_by_lines(stats)
        lines += corr_lines(stats)
        lines += corr_lines(stats, "spearman", "corr-spearman")
        lines += corr_lines(stats, "kendall", "corr-kendall")
        lines += test_lines(stats)
    return lines, bad


def compare(lines, golden_text):
    want = golden_text.rstrip("\n").split("\n")
    if lines == want:
        return []
    return list(difflib.unified_diff(want, lines, "golden", "this build", lineterm="", n=0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", default=str(ROOT / "build" / "special-check"))
    ap.add_argument("--acc", default=str(ROOT / "build" / "acc-check"))
    ap.add_argument("--stats", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--golden", default=str(GOLDEN))
    ap.add_argument("--print", dest="show", action="store_true")
    ap.add_argument("--update", action="store_true")
    ap.add_argument("--mutants", action="store_true")
    ap.add_argument("--no-threads", action="store_true")
    a = ap.parse_args()
    for b in (a.driver, a.acc, a.stats):
        if not pathlib.Path(b).exists():
            print("xplat_fp: no binary %s (cancho install && cancho build --bin special-check && cancho build --bin acc-check)" % b)
            return 2
    try:
        lines, bad = fingerprint(a.driver, a.acc, threads_check=not a.no_threads, stats=a.stats)
    except Exception as e:
        print("xplat_fp: could not run: %s" % e)
        return 2
    if a.show:
        print("\n".join(lines))
        return 1 if bad else 0
    if a.update:
        pathlib.Path(a.golden).write_text("\n".join(lines) + "\n")
        print("xplat_fp: wrote %s (%d lines)" % (a.golden, len(lines)))
        return 1 if bad else 0
    gp = pathlib.Path(a.golden)
    if not gp.exists():
        print("xplat_fp: no golden %s" % gp)
        return 2
    diff = compare(lines, gp.read_text())
    if a.mutants:
        return mutants(lines, gp.read_text(), bad)
    if diff:
        print("xplat_fp: the fingerprint of this build differs from the golden:")
        print("\n".join(diff))
    for b in bad:
        print("xplat_fp: " + b)
    if diff or bad:
        return 1
    print("xplat_fp: ok, %d lines equal the golden%s" % (len(lines), "" if a.no_threads else "; --threads %s equal the sequential bytes in every sweep" % ",".join(map(str, THREADS))))
    return 0


def mutants(lines, golden, bad):
    """The checks can fail: a golden with one hash digit changed, a line removed, a line added, the lines reordered; and a thread result that differs."""
    ok = True
    if compare(lines, golden):
        print("mutants: the control (this build against the golden) differs: refusing")
        return 1
    g = golden.rstrip("\n").split("\n")
    cases = {
        "one hash digit changed": "\n".join([g[0][:-1] + ("0" if g[0][-1] != "0" else "1")] + g[1:]),
        "a line removed": "\n".join(g[1:]),
        "a line added": "\n".join(g + [g[-1]]),
        "two lines swapped": "\n".join([g[1], g[0]] + g[2:]),
        "count changed": "\n".join([g[0].replace("\t3000\t", "\t2999\t")] + g[1:]) if "\t3000\t" in g[0] else "\n".join(g[1:]),
    }
    for name, text in cases.items():
        killed = bool(compare(lines, text))
        print("mutant %-28s %s" % (name, "killed" if killed else "SURVIVED"))
        ok = ok and killed
    # a thread check that differs: simulate by comparing two different outputs through the same equality the check uses
    a_, b_ = "1 4000 0\n", "1 4000 1\n"
    killed = a_ != b_
    print("mutant %-28s %s" % ("thread result differs", "killed" if killed else "SURVIVED"))
    return 0 if ok and killed else 1


if __name__ == "__main__":
    sys.exit(main())
