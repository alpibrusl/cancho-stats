#!/usr/bin/env python3
"""S3.5 gate T-G (docs/test.md sections 4, 5 and 12) for `--correct`:

    python3 scripts/gate_test5.py [--bin build/stats] [--files 100] [--seed 1]

Families of 1 to 60 entries made by every kind of test (`t1` on many columns with constants and equal columns, so that there are null p-values, ties and p-values of 1; `paired`, `t2`, `chisq` and `prop` on several columns),
run with each of `--correct bonferroni`, `holm` and `bh`: `family.p_adjusted` against scripts/oracle_test.py `adjust`, which works from `Fraction` of the printed p-values (BIT FOR BIT), with `size`, `correction` and the
header's `adjustment`; `p_value` of every entry unchanged by the correction (the answer with `--correct none` is the same bytes but for the header and the family). Properties of the program's own answers: p <= p~ <= 1,
the order of the p-values kept, Holm <= Bonferroni, BH <= Holm, a family of one changes nothing, null for an entry without a p-value. CONTROLS (variants of the definitions) must differ from the program on a tenth of the
families: Holm with m - j for m - j + 1, BH without the running minimum, Bonferroni without the cap at 1, Holm without the running maximum. (Ties: the order of equal p-values cannot change an adjusted value, the running maximum and minimum make tied entries equal, so no control exists for it; section 12 says so.)

Exit 0: every number is the expected one. Exit 1: a difference. Exit 2: the harness failed.
"""
import argparse
import json
import pathlib
import subprocess
import sys
import tempfile
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402
import gate_describe as gd  # noqa: E402
import oracle_test as ot  # noqa: E402

METHODS = ("bonferroni", "holm", "bh")


def run(binary, tmp, text, argv):
    (tmp / "f.csv").write_text(text)
    return subprocess.run([binary, "test"] + argv + ["--root", str(tmp), "f.csv"], capture_output=True, text=True, timeout=300)


def gen_t1(r, ncols, n):
    """a table of numeric columns: some constant (no p-value), some equal to another (ties), some with the mean equal to mu0 (p = 1), some far from it (tiny p)"""
    cols = []
    for j in range(ncols):
        kind = int(r.next() % 8)
        if kind == 0:
            col = ["3"] * n
        elif kind == 1 and cols:
            col = list(cols[int(r.next() % len(cols))])
        elif kind == 2:
            col = ["%d" % (3 if i % 2 else 3) if False else ("2" if i % 2 else "4") for i in range(n)]   # mean exactly 3
        elif kind == 3:
            col = ["%r" % (100.0 + (r.next() % 1000) / 1000.0) for _ in range(n)]
        else:
            col = ["%r" % (3.0 + ((r.next() % 20001) - 10000) / 4000.0) for _ in range(n)]
        cols.append(col)
    names = ["c%d" % j for j in range(ncols)]
    return names, cols


def p_values(doc):
    return [x.get("p_value") for x in doc["tests"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--files", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    r = SM64(a.seed)
    bad = []
    ctl = {"holm-m-minus-j": [0, 0], "bh-no-running-min": [0, 0], "bonferroni-uncapped": [0, 0], "holm-no-running-max": [0, 0]}
    checked = 0
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for fi in range(a.files):
            ncols = [1, 2, 3, 5, 10, 20, 60][int(r.next() % 7)]
            n = [3, 5, 12][int(r.next() % 3)]
            names, cols = gen_t1(r, ncols, n)
            text = ",".join(names) + "\n" + "\n".join(",".join(c[i] for c in cols) for i in range(n)) + "\n"
            base = ["t1", "--cols", ",".join(names), "--mu0", "3"]
            p0 = run(a.bin, tmp, text, base)
            if p0.returncode != 0:
                print("gate_test5: exit %d on file %d: %s" % (p0.returncode, fi, p0.stdout[:300]))
                return 2
            d0 = json.loads(p0.stdout)["data"]
            ps = p_values(d0)
            adj_by = {}
            for m in METHODS:
                p = run(a.bin, tmp, text, base + ["--correct", m])
                if p.returncode != 0:
                    print("gate_test5: exit %d on file %d with %s: %s" % (p.returncode, fi, m, p.stdout[:300]))
                    return 2
                d = json.loads(p.stdout)["data"]
                checked += 1
                tag = "file %d %s" % (fi, m)
                fam = d.get("family")
                want = ot.adjust(ps, m)
                live = sum(1 for x in ps if x is not None)
                if d.get("adjustment") != m or fam is None or fam.get("size") != live or fam.get("correction") != m:
                    bad.append((fi, "%s header/family: adjustment %r family %r (want size %d)" % (tag, d.get("adjustment"), {k: v for k, v in (fam or {}).items() if k != "p_adjusted"}, live)))
                    continue
                got = fam["p_adjusted"]
                if len(got) != len(ps) or any(not (gd.same(g, w) if w is not None else g is None) for g, w in zip(got, want)):
                    bad.append((fi, "%s p_adjusted: want %r got %r (p %r)" % (tag, want, got, ps)))
                if p_values(d) != ps:
                    bad.append((fi, "%s: p_value changed by the correction" % tag))
                for q, g in zip(ps, got):
                    if q is None:
                        continue
                    if not (q <= g <= 1.0):
                        bad.append((fi, "%s: p %r adjusted %r outside [p, 1]" % (tag, q, g)))
                live_idx = [i for i, q in enumerate(ps) if q is not None]
                for i in live_idx:
                    for j in live_idx:
                        if ps[i] < ps[j] and not got[i] <= got[j]:
                            bad.append((fi, "%s: the order of the p-values is not kept" % tag))
                            break
                if live == 1 and live_idx and got[live_idx[0]] != ps[live_idx[0]]:
                    bad.append((fi, "%s: a family of one changed the p-value" % tag))
                adj_by[m] = got
            if len(adj_by) == 3:
                for i, q in enumerate(ps):
                    if q is not None and not (adj_by["bh"][i] <= adj_by["holm"][i] <= adj_by["bonferroni"][i]):
                        bad.append((fi, "file %d: BH <= Holm <= Bonferroni violated at %d: %r" % (fi, i, (adj_by["bh"][i], adj_by["holm"][i], adj_by["bonferroni"][i]))))
                        break
                # the controls
                live_idx = [i for i, q in enumerate(ps) if q is not None]
                m_ = len(live_idx)
                if m_ >= 2:
                    exact = {i: Fraction(ps[i]) for i in live_idx}
                    order = sorted(live_idx, key=lambda i: (exact[i], i))
                    run_ = Fraction(0)
                    wrong = [None] * len(ps)
                    for rr, i in enumerate(order, start=1):
                        run_ = max(run_, min(Fraction(1), (m_ - rr) * exact[i]))
                        wrong[i] = ot.rounded(run_)
                    ctl["holm-m-minus-j"][1] += 1
                    ctl["holm-m-minus-j"][0] += 0 if wrong == adj_by["holm"] else 1
                    wrong = [None] * len(ps)
                    for rr, i in enumerate(order, start=1):
                        wrong[i] = ot.rounded(min(Fraction(1), m_ * exact[i] / rr))
                    ctl["bh-no-running-min"][1] += 1
                    ctl["bh-no-running-min"][0] += 0 if wrong == adj_by["bh"] else 1
                    wrong = [None if q is None else ot.rounded(m_ * Fraction(q)) for q in ps]
                    ctl["bonferroni-uncapped"][1] += 1
                    ctl["bonferroni-uncapped"][0] += 0 if wrong == adj_by["bonferroni"] else 1
                    wrong = [None] * len(ps)
                    for rr, i in enumerate(order, start=1):
                        wrong[i] = ot.rounded(min(Fraction(1), (m_ - rr + 1) * exact[i]))
                    ctl["holm-no-running-max"][1] += 1
                    ctl["holm-no-running-max"][0] += 0 if wrong == adj_by["holm"] else 1
            # the other kinds: the family is the same machinery over the entries of the answer
            if fi % 5 == 0:
                pairs = ",".join("c%d:c%d" % (j, (j + 1) % ncols) for j in range(ncols) if ncols > 1)
                if pairs:
                    for m in METHODS:
                        p = run(a.bin, tmp, text, ["paired", "--pairs", pairs, "--correct", m])
                        d = json.loads(p.stdout)["data"] if p.returncode == 0 else None
                        if d is None:
                            bad.append((fi, "paired --correct %s: exit %d" % (m, p.returncode)))
                            continue
                        want = ot.adjust(p_values(d), m)
                        got = d["family"]["p_adjusted"]
                        if d.get("adjustment") != m or any(not (gd.same(g, w) if w is not None else g is None) for g, w in zip(got, want)):
                            bad.append((fi, "paired --correct %s: want %r got %r (adjustment %r)" % (m, want, got, d.get("adjustment"))))
                cat = ["x%d" % (r.next() % 3) for _ in range(n * 4)]
                ctext = "a,b,c,g\n" + "\n".join("%s,%s,%s,%s" % (cat[i], cat[(i * 7 + 1) % len(cat)], "y" if r.next() % 2 else "n", "u" if i % 2 else "v") for i in range(n * 4)) + "\n"
                for kind, argv in (("chisq", ["chisq", "--cols", "a,b,c", "--vs", "g"]), ("prop", ["prop", "--cols", "a,b,c", "--success", "x1", "--p0", "0.3"]), ("prop2", ["prop", "--cols", "a,b,c", "--success", "x1", "--group", "g"])):
                    if argv is None:
                        continue
                    for m in METHODS:
                        p = run(a.bin, tmp, ctext, argv + ["--correct", m])
                        d = json.loads(p.stdout)["data"] if p.returncode == 0 else None
                        if d is None:
                            bad.append((fi, "%s --correct %s: exit %d: %s" % (kind, m, p.returncode, p.stdout[:200])))
                            continue
                        want = ot.adjust(p_values(d), m)
                        got = d["family"]["p_adjusted"]
                        if d.get("adjustment") != m or any(not (gd.same(g, w) if w is not None else g is None) for g, w in zip(got, want)):
                            bad.append((fi, "%s --correct %s: want %r got %r (adjustment %r)" % (kind, m, want, got, d.get("adjustment"))))
    for k, (outside, total) in sorted(ctl.items()):
        if total and outside * 10 < total:
            bad.append((-1, "control %s differed from the program on only %d of %d families: the gate cannot tell it from the program" % (k, outside, total)))
    print("gate_test5: %d files, %d corrected answers compared with the oracle, %d differ; controls that differ from the program: %s" % (a.files, checked, len(bad), ", ".join("%s %d of %d" % (k, v[0], v[1]) for k, v in sorted(ctl.items()))))
    for b in bad[:8]:
        print("  file %d: %s" % b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
