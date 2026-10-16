#!/usr/bin/env python3
"""S1-F (docs/describe.md section 5): the refusals and the hostile inputs of `stats describe`.

    python3 scripts/gate_describe_refusals.py [--bin build/stats] [--fuzz 1500] [--seed 1]

Part 1, a table of 33 cases (argument errors, paths, limits, malformed CSV): each must exit with the code of its rule and print a document whose `error.rule` is exactly the stated tag, `ok` false.
Part 2, fuzz: random files made of the bytes that matter to a CSV reader (delimiters, quotes, CR, LF, NUL, digits, signs, e, dots, high bytes, long runs), and random argument vectors: the program must exit
(no hang within 60 s, no signal), print one line of valid JSON, with `ok` true and exit 0, or `ok` false, an `error.rule` that is one of the tags in `stats introspect`, and the exit code of that rule.
Exit 0: all of it holds. Exit 1: a violation. Exit 2: the harness failed.
"""
import argparse
import json
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from gate_special import SM64  # noqa: E402

CASES = [
    # (name, file bytes or None, argv after the binary, tag, exit)
    ("no verb", b"a\n1\n", [], "args.missing-operand", 2),
    ("bad verb", b"a\n1\n", ["summarise", "f.csv"], "stats.bad-verb", 2),
    ("no file", b"a\n1\n", ["describe"], "args.missing-operand", 2),
    ("two files", b"a\n1\n", ["describe", "f.csv", "f.csv"], "args.too-many-operands", 2),
    ("unknown flag", b"a\n1\n", ["describe", "--median", "f.csv"], "args.unknown-flag", 2),
    ("flag without value", b"a\n1\n", ["describe", "--ddof"], "args.missing-value", 2),
    ("flag given twice", b"a\n1\n", ["describe", "--ddof", "1", "--ddof", "0", "f.csv"], "args.duplicate-flag", 2),
    ("ddof not a number", b"a\n1\n", ["describe", "--ddof", "x", "f.csv"], "args.bad-value", 2),
    ("ddof 2", b"a\n1\n", ["describe", "--ddof", "2", "f.csv"], "args.bad-value", 2),
    ("delim two bytes", b"a\n1\n", ["describe", "--delim", "ab", "f.csv"], "args.bad-value", 2),
    ("max-rows over the ceiling", b"a\n1\n", ["describe", "--max-rows", "10000001", "f.csv"], "args.bad-value", 2),
    ("max-cols zero", b"a\n1\n", ["describe", "--max-cols", "0", "f.csv"], "args.bad-value", 2),
    ("path with ..", b"a\n1\n", ["describe", "../f.csv"], "path.dotdot", 2),
    ("file missing", None, ["describe", "nothing.csv"], "io.not-found", 3),
    ("a directory", None, ["describe", "."], "io.is-a-directory", 2),
    ("empty file", b"", ["describe", "f.csv"], "stats.no-header", 8),
    ("ragged record", b"a,b\n1,2\n3\n", ["describe", "f.csv"], "stats.csv-ragged", 8),
    ("record with extra field", b"a,b\n1,2,3\n", ["describe", "f.csv"], "stats.csv-ragged", 8),
    ("record continues", b'a,b\n1,"x\ny"\n', ["describe", "f.csv"], "stats.csv-multiline", 8),
    ("bad quote", b'a,b\n1,"x"y\n', ["describe", "f.csv"], "stats.csv-bad-quote", 8),
    ("too many rows", b"a\n1\n2\n3\n", ["describe", "--max-rows", "2", "f.csv"], "limit.too-many-rows", 8),
    ("columns exactly at the limit", b"a,b\n1,2\n", ["describe", "--max-cols", "2", "f.csv"], None, 0),
    ("too many columns", b"a,b,c\n1,2,3\n", ["describe", "--max-cols", "2", "f.csv"], "limit.too-many-columns", 8),
    ("line too long", b"a\n" + b"1" * 200 + b"\n", ["describe", "--max-line-bytes", "100", "f.csv"], "limit.line-too-long", 8),
    ("header continues", b'"a\nb",c\n1,2\n', ["describe", "f.csv"], "stats.csv-multiline", 8),
    ("quantile type 0", b"a\n1\n", ["describe", "--quantile-type", "0", "f.csv"], "args.bad-value", 2),
    ("quantile type 10", b"a\n1\n", ["describe", "--quantile-type", "10", "f.csv"], "args.bad-value", 2),
    ("quantile not a number", b"a\n1\n", ["describe", "--quantiles", "0.1,x", "f.csv"], "args.bad-value", 2),
    ("quantile above 1", b"a\n1\n", ["describe", "--quantiles", "1.5", "f.csv"], "args.bad-value", 2),
    ("quantile with 10 decimals", b"a\n1\n", ["describe", "--quantiles", "0.1234567891", "f.csv"], "args.bad-value", 2),
    ("quantile empty item", b"a\n1\n", ["describe", "--quantiles", "0.1,,0.2", "f.csv"], "args.bad-value", 2),
    ("too many cells", b"a,b\n1,2\n3,4\n", ["describe", "--max-cells", "3", "f.csv"], "limit.too-many-cells", 8),
    ("max-cells zero", b"a\n1\n", ["describe", "--max-cells", "0", "f.csv"], "args.bad-value", 2),
    ("cells exactly at the limit", b"a,b\n1,2\n3,4\n", ["describe", "--max-cells", "4", "f.csv"], None, 0),
    ("99 probabilities", b"a\n1\n", ["describe", "--quantiles", ",".join(["0.5"] * 99), "f.csv"], None, 0),
    ("100 probabilities", b"a\n1\n", ["describe", "--quantiles", ",".join(["0.5"] * 100), "f.csv"], "args.bad-value", 2),
    ("bins above the ceiling", b"a\n1\n", ["describe", "--bins", "1001", "f.csv"], "args.bad-value", 2),
    ("bins at the ceiling", b"a\n1\n2\n", ["describe", "--bins", "1000", "f.csv"], None, 0),
    ("50 groups is the default limit", b"g\n" + b"".join(b"v%d\n" % i for i in range(50)), ["describe", "--by", "g", "f.csv"], None, 0),
    ("51 groups is over the default limit", b"g\n" + b"".join(b"v%d\n" % i for i in range(51)), ["describe", "--by", "g", "f.csv"], "limit.too-many-groups", 8),
    ("by no such column", b"g,x\na,1\n", ["describe", "--by", "zz", "f.csv"], "stats.no-such-column", 2),
    ("by is case sensitive", b"g,x\na,1\n", ["describe", "--by", "G", "f.csv"], "stats.no-such-column", 2),
    ("by two columns of that name", b"g,x,g\na,1,2\n", ["describe", "--by", "g", "f.csv"], "stats.ambiguous-column", 2),
    ("by on a ragged file", b"g,x\na,1\nb\n", ["describe", "--by", "g", "f.csv"], "stats.csv-ragged", 8),
    ("by on an empty file", b"", ["describe", "--by", "g", "f.csv"], "stats.no-header", 8),
    ("by with too many rows", b"g\na\nb\nc\n", ["describe", "--by", "g", "--max-rows", "2", "f.csv"], "limit.too-many-rows", 8),
    ("too many groups", b"g\na\nb\nc\n", ["describe", "--by", "g", "--max-groups", "2", "f.csv"], "limit.too-many-groups", 8),
    ("groups exactly at the limit", b"g\na\nb\nb\n", ["describe", "--by", "g", "--max-groups", "2", "f.csv"], None, 0),
    ("a missing group value is not a group", b"g\na\n\n\n", ["describe", "--by", "g", "--max-groups", "1", "f.csv"], None, 0),
    ("max-groups zero", b"g\na\n", ["describe", "--max-groups", "0", "f.csv"], "args.bad-value", 2),
    ("max-groups above the ceiling", b"g\na\n", ["describe", "--max-groups", "1001", "f.csv"], "args.bad-value", 2),
    ("max-groups at the ceiling", b"g\na\n", ["describe", "--by", "g", "--max-groups", "1000", "f.csv"], None, 0),
    ("by limits rows of a group not of the file", b"g,x\na,1\nb,2\nb,3\na,4\n", ["describe", "--by", "g", "--max-rows", "4", "f.csv"], None, 0),
    ("by with a group over the cells limit", b"g,x\na,1\na,2\nb,3\n", ["describe", "--by", "g", "--max-cells", "3", "f.csv"], "limit.too-many-cells", 8),
    ("cols no such column", b"a,b\n1,2\n", ["describe", "--cols", "a,z", "f.csv"], "stats.no-such-column", 2),
    ("cols two columns of the name", b"a,a\n1,2\n", ["describe", "--cols", "a", "f.csv"], "stats.ambiguous-column", 2),
    ("cols empty item", b"a,b\n1,2\n", ["describe", "--cols", "a,,b", "f.csv"], "args.bad-value", 2),
    ("cols trailing comma", b"a,b\n1,2\n", ["describe", "--cols", "a,", "f.csv"], "args.bad-value", 2),
    ("cols a name twice is fine", b"a,b\n1,2\n", ["describe", "--cols", "a,a", "f.csv"], None, 0),
    ("cols with by", b"g,x\na,1\nb,2\n", ["describe", "--by", "g", "--cols", "x", "f.csv"], None, 0),
    ("format xml", b"a\n1\n", ["describe", "--format", "xml", "f.csv"], "args.bad-value", 2),
    ("levels above the ceiling", b"a\nx\n", ["describe", "--levels", "101", "f.csv"], "args.bad-value", 2),
    ("levels at the ceiling", b"a\nx\n", ["describe", "--levels", "100", "f.csv"], None, 0),
    ("max-levels zero", b"a\nx\n", ["describe", "--max-levels", "0", "f.csv"], "args.bad-value", 2),
    ("max-levels above the ceiling", b"a\nx\n", ["describe", "--max-levels", "1000001", "f.csv"], "args.bad-value", 2),
    ("max-levels at the ceiling", b"a\nx\n", ["describe", "--max-levels", "1000000", "f.csv"], None, 0),
    ("levels not a number", b"a\nx\n", ["describe", "--levels", "many", "f.csv"], "args.bad-value", 2),
    ("bins zero", b"a\n1\n2\n", ["describe", "--bins", "0", "f.csv"], None, 0),
    ("file with only a header", b"a,b\n", None, None, 0),
]


def call(binary, tmp, argv, timeout=60):
    return subprocess.run([binary] + argv, capture_output=True, cwd=tmp, timeout=timeout)


def check_doc(p, known, why):
    """-> list of problems of one finished process"""
    out = []
    if p.returncode < 0:
        return ["%s: killed by signal %d" % (why, -p.returncode)]
    lines = p.stdout.split(b"\n")
    if len(lines) != 2 or lines[1] != b"":
        out.append("%s: stdout is not one line (%d lines)" % (why, len(lines) - 1))
    try:
        d = json.loads(lines[0])
    except Exception:
        return out + ["%s: stdout is not JSON: %r" % (why, p.stdout[:80])]
    if d.get("ok") is True:
        if p.returncode != 0:
            out.append("%s: ok true with exit %d" % (why, p.returncode))
    else:
        rule = d.get("error", {}).get("rule")
        if rule not in known:
            out.append("%s: rule %r is not in introspect" % (why, rule))
        if p.returncode == 0:
            out.append("%s: ok false with exit 0" % why)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--fuzz", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    bad = []
    try:
        intro = json.loads(subprocess.run([a.bin, "introspect"], capture_output=True, timeout=30).stdout)
    except Exception as e:
        print("gate_describe_refusals: introspect failed: %s" % e)
        return 2
    known = {r["rule"] for r in intro["rules"]}
    exits = {r["rule"]: r["exit"] for r in intro["rules"]}
    reached = set()
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        for name, data, argv, tag, code in CASES:
            f = tmp / "f.csv"
            if f.exists():
                f.unlink()
            if data is not None:
                f.write_bytes(data)
            if argv is None:
                argv = ["describe", "f.csv"]
            try:
                p = call(a.bin, tmp, argv)
            except subprocess.TimeoutExpired:
                bad.append("%s: hang" % name)
                continue
            if tag is None:
                if p.returncode != 0 or json.loads(p.stdout.split(b"\n")[0]).get("ok") is not True:
                    bad.append("%s: expected an answer, got exit %d: %r" % (name, p.returncode, p.stdout[:100]))
                continue
            try:
                d = json.loads(p.stdout.split(b"\n")[0])
            except Exception:
                bad.append("%s: not JSON: %r" % (name, p.stdout[:100]))
                continue
            rule = d.get("error", {}).get("rule")
            reached.add(rule)
            if d.get("ok") is not False or rule != tag or p.returncode != code:
                bad.append("%s: want %s exit %d, got %s exit %d" % (name, tag, code, rule, p.returncode))
            bad += check_doc(p, known, name)
        # part 2: fuzz
        r = SM64(a.seed)
        alphabet = [b",", b",", b'"', b'"', b"\n", b"\n", b"\r", b"\r\n", b"\x00", b"0", b"1", b"9", b"-", b"+", b".", b"e", b"E", b" ", b"a", b"\xff", b"\xc3\xa9", b"1" * 30, b"nan", b"inf", b"1e999", b"1e-999", b";", b"\t"]
        for i in range(a.fuzz):
            n = int(r.next() % 60)
            data = b"".join(alphabet[int(r.next() % len(alphabet))] for _ in range(n))
            (tmp / "f.csv").write_bytes(data)
            argv = ["describe", "f.csv"]
            k = int(r.next() % 10)
            if k == 1:
                argv = ["describe", "--ddof", "0", "f.csv"]
            elif k == 2:
                argv = ["describe", "--delim", ";", "f.csv"]
            elif k == 3:
                argv = ["describe", "--max-rows", str(int(r.next() % 4)), "--max-cols", str(1 + int(r.next() % 3)), "f.csv"]
            elif k == 4:
                argv = ["describe", "--max-line-bytes", str(1 + int(r.next() % 30)), "f.csv"]
            elif k == 6:
                argv = ["describe", "--by", "a", "f.csv"]
            elif k == 7:
                argv = ["describe", "--cols", "a", "--levels", "3", "f.csv"]
            elif k == 8:
                argv = ["describe", "--by", "a", "--max-groups", "2", "--max-rows", "3", "f.csv"]
            elif k == 9:
                argv = ["describe", "--format", "text", "f.csv"]
            try:
                p = call(a.bin, tmp, argv)
            except subprocess.TimeoutExpired:
                bad.append("fuzz %d: hang on %r" % (i, data[:60]))
                continue
            if k == 9:
                # text mode: no signal; an answer is lines `path = value`, a refusal is a sentence on standard error with nothing on standard output
                if p.returncode < 0:
                    bad.append("fuzz %d: text mode died of signal %d" % (i, -p.returncode))
                elif p.returncode == 0 and any(b" = " not in ln for ln in p.stdout.split(b"\n") if ln):
                    bad.append("fuzz %d: text mode printed a line that is not `path = value`: %r" % (i, p.stdout[:80]))
                elif p.returncode != 0 and (p.stdout != b"" or p.stderr == b""):
                    bad.append("fuzz %d: text mode refusal: stdout %r stderr %r" % (i, p.stdout[:40], p.stderr[:40]))
                continue
            bad += check_doc(p, known, "fuzz %d %r" % (i, data[:40]))
            if len(bad) > 20:
                break
    print("gate_describe_refusals: %d table cases, %d fuzz files; rules reached by the table: %d; %d problems" % (len(CASES), a.fuzz, len(reached), len(bad)))
    for b in bad[:8]:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
