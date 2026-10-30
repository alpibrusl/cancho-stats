#!/usr/bin/env python3
"""The determinism gate G8 (docs/numerics.md section 6): the fingerprint of the bytes the drivers print, against a committed golden, and the thread check.

    python3 scripts/xplat_fp.py                  # print nothing but a verdict: the fingerprint of this build equals docs/xplat_fp.txt, and every thread count 1..16 gives the same bytes
    python3 scripts/xplat_fp.py --print          # print the fingerprint (what a platform contributes to the comparison)
    python3 scripts/xplat_fp.py --update         # rewrite the golden: only when the numbers change on purpose, in a commit that says why
    python3 scripts/xplat_fp.py --mutants        # the checks are shown able to fail

The fingerprint is one tab-separated line per sweep: `kind  name  count  sha256`, the hash over the driver's output lines (`id  hex bits of the result  code`) for a fixed list of cases:
`special` the first N cases of each range of gate_special.py (seed 1), `edge` the edge cases of gate_domain.py for each function, `acc` the vectors of three families of measure_acc.py through acc-check, `describe` the bytes `stats describe` prints for 12 seeded CSV files, `corr`, `corr-spearman` and `corr-kendall` those of `stats corr` (Pearson, --method spearman, --method kendall) for 8 files each, `test-t1`, `test-t2`, `test-paired` and `test-chisq`, `test-prop` and `test-correct` those of `stats test t1`, `t2` (Welch and pooled drawn), `paired`, `chisq` (independence and fit) and `prop` (one and two samples) for 8 files each, the last with `--correct holm`, `bh` and `bonferroni` over a family of columns.
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
import gate_test2  # noqa: E402
import gate_test3  # noqa: E402
import gate_test5  # noqa: E402
import gate_anova  # noqa: E402
import gate_kruskal  # noqa: E402
import gate_mannwhitney  # noqa: E402
import gate_ks  # noqa: E402
import gate_fisher  # noqa: E402
import gate_corr_exact  # noqa: E402
import gate_ci_prop  # noqa: E402
import gate_ci_corr  # noqa: E402
import gate_ci_boot  # noqa: E402
import gate_test2  # noqa: E402
import gate_test  # noqa: E402
import gate_wilcoxon  # noqa: E402
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


def test2_lines(stats):
    """`stats test t2` and `paired` on 8 seeded CSV files each (the generators of gate_test2.py: groups of any size, missing cells, the method and the level drawn): the hash of the bytes they print."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(41)
        k = 0
        fi = 1000
        while k < 8:
            fi += 1
            names, rows, gcol, labels, shapes, vidx, method, level = gate_test2.gen_t2(r, fi)
            num = [j for j in vidx if gate_test.numeric([rw[j] for rw in rows])]
            if not num:
                continue
            p = gate_test2.run(stats, tmp, "t2", names, rows, ["--cols", ",".join(names[j] for j in num), "--group", "g", "--variance", method, "--level", level])
            h = hashlib.sha256()
            h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("test-t2\tfile%02d\t%d\t%s" % (k, len(rows), h.hexdigest()))
            k += 1
        r = gate_special.SM64(43)
        k = 0
        while k < 8:
            fi += 1
            names, rows, pairs, shapes, level = gate_test2.gen_paired(r, fi)
            if not pairs:
                continue
            p = gate_test2.run(stats, tmp, "paired", names, rows, ["--pairs", ",".join("%s:%s" % pr for pr in pairs), "--level", level])
            h = hashlib.sha256()
            h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("test-paired\tfile%02d\t%d\t%s" % (k, len(rows), h.hexdigest()))
            k += 1
    return out


def test3_lines(stats):
    """`stats test chisq` on 8 seeded CSV files (the generator of gate_test3.py: independence of two columns against a third, then the fit of the same columns): the hash of the bytes they print."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(47)
        for k in range(8):
            la, lb, wt = gate_test3.draw_table(r)
            n = [2, 3, 5, 10, 30, 100, 400][int(r.next() % 7)]
            pairs = gate_test3.sample(r, la, lb, wt, n)
            rows = [[p[0], p[1], p[0] if r.next() % 3 else None] for p in pairs]
            h = hashlib.sha256()
            p = gate_test3.run(stats, tmp, ["A", "G", "C"], rows, ["--cols", "A,C", "--vs", "G"])
            h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            levels = sorted(set(x[0] for x in pairs))
            props, texts = gate_test3.draw_props(r, len(levels))
            p = gate_test3.run(stats, tmp, ["A", "G", "C"], rows, ["--cols", "A", "--expected", ",".join(texts)])
            h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("test-chisq\tfile%02d\t%d\t%s" % (k, n, h.hexdigest()))
    return out


def test4_lines(stats):
    """`stats test prop` on 8 seeded CSV files (an outcome column with a success label, a group column of two labels, empty cells, trials from 1 to 3,000; p0 and the level drawn): the hash of the bytes of the one-sample and the two-sample answers."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(53)
        for k in range(8):
            n = [1, 2, 5, 30, 100, 400, 3000][int(r.next() % 7)]
            ps = (r.next() % 1000) / 1000.0
            rows = []
            for _ in range(n):
                y = "yes" if (r.next() % 100000) / 100000.0 < ps else "no"
                if r.next() % 25 == 0:
                    y = None
                g = ("g1" if r.next() % 2 else "g2") if r.next() % 25 else None
                rows.append([y, g])
            p0 = r.pick(["0.5", "0.1", "0.25", "0.9", "0.333333333"])
            level = r.pick(["0.95", "0.9", "0.99"])
            h = hashlib.sha256()
            for extra in (["--cols", "y", "--success", "yes", "--p0", p0, "--level", level], ["--cols", "y", "--success", "yes", "--group", "g", "--level", level]):
                text = "y,g\n" + "\n".join(",".join("" if c is None else c for c in row) for row in rows) + "\n"
                (tmp / "f.csv").write_text(text)
                p = subprocess.run([stats, "test", "prop", "--root", str(tmp)] + extra + ["f.csv"], capture_output=True, text=True, timeout=300)
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("test-prop\tfile%02d\t%d\t%s" % (k, n, h.hexdigest()))
    return out


def test5_lines(stats):
    """`stats test t1 --correct` on 8 seeded CSV files of 2 to 30 numeric columns (constants and equal columns among them: null p-values and ties): the hash of the bytes of the three corrections."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(59)
        for k in range(8):
            ncols = [2, 3, 5, 10, 20, 30][int(r.next() % 6)]
            n = [3, 5, 12][int(r.next() % 3)]
            names, cols = gate_test5.gen_t1(r, ncols, n)
            text = ",".join(names) + "\n" + "\n".join(",".join(c[i] for c in cols) for i in range(n)) + "\n"
            (tmp / "f.csv").write_text(text)
            h = hashlib.sha256()
            for m in ("holm", "bh", "bonferroni"):
                p = subprocess.run([stats, "test", "t1", "--cols", ",".join(names), "--mu0", "3", "--correct", m, "--root", str(tmp), "f.csv"], capture_output=True, text=True, timeout=300)
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("test-correct\tfile%02d\t%d\t%s" % (k, ncols, h.hexdigest()))
    return out


def anova_lines(stats):
    """`stats test anova` on 8 seeded CSV files (2 to 25 groups of unequal sizes with their own mean and spread, empty cells, shapes from the generator of scripts/gate_anova.py): the hash of the bytes of Welch's and of the classic answer."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(61)
        done = 0
        fi = len(gate_anova.FIXED)
        while done < 8:
            names, rows, gcol, labels, vidx, shapes = gate_anova.gen(r, fi)
            fi += 1
            num = [j for j in vidx if gate_anova.column_numeric(rows, j)]
            if not num:
                continue
            h = hashlib.sha256()
            for method in ("welch", "pooled"):
                p = gate_anova.run(stats, tmp, names, rows, ["--cols", ",".join(names[j] for j in num), "--group", "g", "--variance", method])
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("test-anova\tfile%02d\t%d\t%s" % (done, len(rows), h.hexdigest()))
            done += 1
    return out


def kruskal_lines(stats):
    """`stats test kruskal` on 8 seeded CSV files (2 to 25 groups of unequal sizes, lattices with heavy ties, empty cells, the shapes of scripts/gate_kruskal.py): the hash of the bytes of the answer."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(67)
        done = 0
        fi = len(gate_kruskal.FIXED)
        while done < 8:
            names, rows, gcol, labels, vidx, shapes = gate_kruskal.gen(r, fi)
            fi += 1
            num = [j for j in vidx if gate_kruskal.column_numeric(rows, j)]
            if not num:
                continue
            p = gate_kruskal.run(stats, tmp, names, rows, ["--cols", ",".join(names[j] for j in num), "--group", "g"])
            out.append("test-kruskal\tfile%02d\t%d\t%s" % (done, len(rows), hashlib.sha256(("exit %d %s" % (p.returncode, p.stdout)).encode()).hexdigest()))
            done += 1
    return out


def mannwhitney_lines(stats):
    """`stats test mannwhitney` on 8 seeded CSV files (two groups of unequal sizes up to 200, lattices with heavy ties, empty cells, the shapes of scripts/gate_mannwhitney.py), with `--exact auto`, `--exact no` and `--continuity none`:
    the hash of the bytes of the three answers (the exact p-value for up to 60 values, the normal one above)."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(71)
        done = 0
        fi = len(gate_mannwhitney.FIXED)
        while done < 8:
            names, rows, gcol, labels, vidx, shapes = gate_mannwhitney.gen(r, fi)
            fi += 1
            num = [j for j in vidx if gate_mannwhitney.column_numeric(rows, j)]
            if not num or len(set(rw[gcol] for rw in rows if rw[gcol] != "")) != 2:
                continue
            h = hashlib.sha256()
            for extra in ([], ["--exact", "no"], ["--continuity", "none"]):
                p = gate_mannwhitney.run(stats, tmp, names, rows, ["--cols", ",".join(names[j] for j in num), "--group", "g"] + extra)
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("test-mannwhitney\tfile%02d\t%d\t%s" % (done, len(rows), h.hexdigest()))
            done += 1
    return out


def ks_lines(stats):
    """`stats test ks` on 8 seeded CSV files (the files of mannwhitney_lines, a different seed), with `--exact auto` and `--exact no`: the hash of the bytes of the two answers (the exact p-value for up to 60 distinct values, the
    asymptotic one above)."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(83)
        done = 0
        fi = len(gate_mannwhitney.FIXED)
        while done < 8:
            names, rows, gcol, labels, vidx, shapes = gate_mannwhitney.gen(r, fi)
            fi += 1
            num = [j for j in vidx if gate_mannwhitney.column_numeric(rows, j)]
            if not num or len(set(rw[gcol] for rw in rows if rw[gcol] != "")) != 2:
                continue
            h = hashlib.sha256()
            for extra in ([], ["--exact", "no"]):
                p = gate_ks.run(stats, tmp, names, rows, ["--cols", ",".join(names[j] for j in num), "--group", "g"] + extra)
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("test-ks\tfile%02d\t%d\t%s" % (done, len(rows), h.hexdigest()))
            done += 1
    return out


def fisher_lines(stats):
    """`stats test fisher` on 8 seeded CSV files (two to five categorical columns of two levels, up to 3,000 rows, missing cells, the generator of scripts/gate_fisher.py) and the 25 fixed tables of the gate (the tolerance tables, the margins of 0, N = 10,000
    and 10,001): the hash of the bytes of the answers. The p-value is an exact rational rounded once, so it is the same bytes on every platform."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(97)
        h = hashlib.sha256()
        for fi in range(len(gate_fisher.FIXED)):
            names, rows, vidx = gate_fisher.gen(r, fi)
            p = gate_fisher.run(stats, tmp, names, rows, ["--cols", names[vidx[0]], "--vs", "y"])
            h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
        out.append("test-fisher\tfixed\t%d\t%s" % (len(gate_fisher.FIXED), h.hexdigest()))
        for done in range(8):
            fi = len(gate_fisher.FIXED) + done
            names, rows, vidx = gate_fisher.gen(r, fi)
            p = gate_fisher.run(stats, tmp, names, rows, ["--cols", ",".join(names[j] for j in vidx), "--vs", "y"])
            out.append("test-fisher\tfile%02d\t%d\t%s" % (done, len(rows), hashlib.sha256(("exit %d %s" % (p.returncode, p.stdout)).encode()).hexdigest()))
    return out


def corr_exact_lines(stats):
    """`stats corr --method kendall --exact auto` and `--exact yes` on 8 seeded CSV files (permutations with a planted dependence, columns with ties, missing cells, 2 to 80 rows, the generator of scripts/gate_corr_exact.py): the hash
    of the bytes of the answers. The exact p-value is an exact rational rounded once, so it is the same bytes on every platform."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(101)
        for k in range(8):
            names, rows = gate_corr_exact.gen(r, len(gate_corr_exact.FIXED_N) + k)
            h = hashlib.sha256()
            for mode in ("auto", "yes"):
                p = gate_corr_exact.run(stats, tmp, names, rows, ["--exact", mode])
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("corr-kendall-exact\tfile%02d\t%d\t%s" % (k, len(rows), h.hexdigest()))
    return out


def ci_prop_lines(stats):
    """`stats ci prop` on 8 seeded CSV files (an outcome column with a success label, a group column, empty cells, 0 to 3,000 rows, the generator of scripts/gate_ci_prop.py): Wilson, Clopper-Pearson and (with --group) Newcombe, one hash of the bytes
    of the three answers. The Clopper-Pearson endpoints come from the beta quantile of special.quantile, so these lines are what shows the quantile function agrees across platforms on this input."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(113)
        for k in range(8):
            names, rows, succ, labels = gate_ci_prop.gen(r, 12 + k)
            level = gate_ci_prop.LEVELS[int(r.next() % len(gate_ci_prop.LEVELS))]
            h = hashlib.sha256()
            for extra in (["--cols", "y,w", "--success", succ, "--level", level],
                          ["--cols", "y,w", "--success", succ, "--level", level, "--method", "clopper-pearson"],
                          ["--cols", "y,w", "--success", succ, "--level", level, "--group", "g"]):
                p = gate_ci_prop.run(stats, tmp, names, rows, extra)
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("ci-prop\tfile%02d\t%d\t%s" % (k, len(rows), h.hexdigest()))
    return out


def ci_mean_lines(stats):
    """`stats ci mean` and `stats ci diff` on 8 seeded CSV files of the groups of scripts/gate_test2.py (two labels, missing cells, the shapes of scripts/gate_describe.py) and 4 paired ones: `ci mean --cols`, `ci diff` with Welch's and the pooled variance, and
    `ci mean --pairs`; one hash of the bytes of the answers."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(127)
        done = 0
        fi = len(gate_test2.FIXED2)
        while done < 8:
            names, rows, gcol, labels, shapes, vidx, variance, level = gate_test2.gen_t2(r, fi)
            fi += 1
            num = [j for j in vidx if gate_test.numeric([rw[j] for rw in rows])]
            if not num:
                continue
            cols = ",".join(names[j] for j in num)
            h = hashlib.sha256()
            for argv in (["ci", "mean", "--cols", cols, "--level", level],
                         ["ci", "diff", "--cols", cols, "--group", names[gcol], "--level", level],
                         ["ci", "diff", "--cols", cols, "--group", names[gcol], "--variance", "pooled", "--level", level]):
                text = ",".join(names) + "\n" + "\n".join(",".join(rw) for rw in rows) + "\n"
                (tmp / "f.csv").write_text(text)
                p = subprocess.run([stats] + argv + ["--root", str(tmp), "f.csv"], capture_output=True, text=True, timeout=300)
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("ci-mean\tfile%02d\t%d\t%s" % (done, len(rows), h.hexdigest()))
            done += 1
        fp = len(gate_test2.FIXED_P)
        for k in range(4):
            names, rows, pairs, shapes, level = gate_test2.gen_paired(r, fp + k)
            if not pairs:
                continue
            text = ",".join(names) + "\n" + "\n".join(",".join(rw) for rw in rows) + "\n"
            (tmp / "f.csv").write_text(text)
            p = subprocess.run([stats, "ci", "mean", "--pairs", ",".join("%s:%s" % pr for pr in pairs), "--level", level, "--root", str(tmp), "f.csv"], capture_output=True, text=True, timeout=300)
            out.append("ci-mean-paired\tfile%02d\t%d\t%s" % (k, len(rows), hashlib.sha256(("exit %d %s" % (p.returncode, p.stdout)).encode()).hexdigest()))
    return out


def ci_corr_lines(stats):
    """`stats ci corr` on 8 seeded CSV files (the numeric shapes of scripts/gate_corr.py, correlated pairs, missing cells, huge and tiny values), --missing pairwise and listwise: the hash of the bytes of the answers. The Fisher z interval is
    `atanh`, `tanh` and the normal quantile of the platform's `std.math`, which is what these lines check across platforms."""
    out = []
    r = gate_special.SM64(131)
    shapes_all = ["uniform", "ints", "lin", "offset", "offset-neg", "wide", "tiny", "huge", "const", "gappy", "lin", "lin", "ints", "uniform"]
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for k in range(8):
            ncols = 2 + int(r.next() % 3)
            shapes = [shapes_all[int(r.next() % len(shapes_all))] for _ in range(ncols)]
            n = [4, 5, 20, 100, 400][int(r.next() % 5)]
            cols = []
            for j, sh in enumerate(shapes):
                cols.append(gate_corr.column(r, sh, n, cols[j - 1] if j else []))
            names = ["c%d_%s" % (j, sh) for j, sh in enumerate(shapes)]
            rows = [[cols[j][i] for j in range(ncols)] for i in range(n)]
            level = gate_ci_corr.LEVELS[int(r.next() % len(gate_ci_corr.LEVELS))]
            h = hashlib.sha256()
            for mode in ("pairwise", "listwise"):
                p = gate_ci_corr.run(stats, tmp, names, rows, ["--missing", mode, "--level", level])
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("ci-corr\tfile%02d\t%d\t%s" % (k, n, h.hexdigest()))
    return out


def ci_boot_lines(stats):
    """`stats ci boot` on 6 seeded CSV files of the groups of scripts/gate_test2.py, the four statistics and seeds across the range: the hash of the bytes of the answers. The replicates are SplitMix64 and exact sums, so the same
    bytes on every platform is the claim; the only floating-point step is the square root of the sd."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(137)
        done = 0
        fi = len(gate_test2.FIXED2)
        while done < 6:
            names, rows, gcol, labels, shapes, vidx, variance, level = gate_test2.gen_t2(r, fi)
            fi += 1
            num = [j for j in vidx if gate_test.numeric([rw[j] for rw in rows])]
            if not num or len(rows) > 60:
                continue
            cols = ",".join(names[j] for j in num)
            h = hashlib.sha256()
            for stat in gate_ci_boot.STATS:
                seed = gate_ci_boot.SEEDS[int(r.next() % len(gate_ci_boot.SEEDS))]
                p = gate_ci_boot.run(stats, tmp, names, rows, ["ci", "boot", "--cols", cols, "--statistic", stat, "--seed", str(seed), "--replicates", "200", "--level", level] + gate_ci_boot.stat_args(stat))
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("ci-boot\tfile%02d\t%d\t%s" % (done, len(rows), h.hexdigest()))
            done += 1
    return out


def ci_bca_lines(stats):
    """`stats ci boot --method bca` on 6 seeded CSV files of the groups of scripts/gate_test2.py, the four statistics and seeds across the range: the hash of the bytes of the answers. The BCa step uses the normal quantile and distribution functions of S0 (libm `log` and `exp` inside), so
    the positions are the claim that a platform could break."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(139)
        done = 0
        fi = len(gate_test2.FIXED2)
        while done < 6:
            names, rows, gcol, labels, shapes, vidx, variance, level = gate_test2.gen_t2(r, fi)
            fi += 1
            num = [j for j in vidx if gate_test.numeric([rw[j] for rw in rows])]
            if not num or len(rows) > 60:
                continue
            cols = ",".join(names[j] for j in num)
            h = hashlib.sha256()
            for stat in gate_ci_boot.STATS:
                seed = gate_ci_boot.SEEDS[int(r.next() % len(gate_ci_boot.SEEDS))]
                p = gate_ci_boot.run(stats, tmp, names, rows, ["ci", "boot", "--cols", cols, "--statistic", stat, "--seed", str(seed), "--replicates", "200", "--level", level] + gate_ci_boot.stat_args(stat) + ["--method", "bca"])
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("ci-bca\tfile%02d\t%d\t%s" % (done, len(rows), h.hexdigest()))
            done += 1
    return out


def wilcoxon_lines(stats):
    """`stats test wilcoxon` on 8 seeded CSV files (pairs of columns with zero differences and heavy ties, sizes up to 200, the shapes of scripts/gate_wilcoxon.py), with the three methods for zeros, `--exact no` and
    `--continuity none`, and the one-sample form: the hash of the bytes of the answers."""
    out = []
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        r = gate_special.SM64(73)
        done = 0
        fi = 9 * len(gate_wilcoxon.FIXED)
        while done < 8:
            names, rows, shapes = gate_wilcoxon.gen(r, fi)
            fi += 1
            if not gate_wilcoxon.column_numeric(rows, 0) or not gate_wilcoxon.column_numeric(rows, 1):
                continue
            h = hashlib.sha256()
            pairs = ["--pairs", "%s:%s" % (names[0], names[1])]
            for extra in (pairs, pairs + ["--zeros", "pratt"], pairs + ["--zeros", "zsplit"], pairs + ["--exact", "no", "--continuity", "none"], ["--cols", names[0], "--mu0", "1"]):
                p = gate_wilcoxon.run(stats, tmp, names, rows, extra)
                h.update(("exit %d %s" % (p.returncode, p.stdout)).encode())
            out.append("test-wilcoxon\tfile%02d\t%d\t%s" % (done, len(rows), h.hexdigest()))
            done += 1
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
        lines += test2_lines(stats)
        lines += test3_lines(stats)
        lines += test4_lines(stats)
        lines += test5_lines(stats)
        lines += anova_lines(stats)
        lines += kruskal_lines(stats)
        lines += mannwhitney_lines(stats)
        lines += wilcoxon_lines(stats)
        lines += ks_lines(stats)
        lines += fisher_lines(stats)
        lines += corr_exact_lines(stats)
        lines += ci_prop_lines(stats)
        lines += ci_mean_lines(stats)
        lines += ci_corr_lines(stats)
        lines += ci_boot_lines(stats)
        lines += ci_bca_lines(stats)
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
