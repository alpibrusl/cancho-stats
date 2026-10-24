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
    # --by reads the header in a pass of its own (discover), with its own copy of the column limit
    ("by: columns exactly at the limit", b"a,b\n1,2\n", ["describe", "--by", "a", "--max-cols", "2", "f.csv"], None, 0),
    ("by: too many columns", b"a,b,c\n1,2,3\n", ["describe", "--by", "a", "--max-cols", "2", "f.csv"], "limit.too-many-columns", 8),
    # corr (docs/corr.md section 4)
    ("corr: two numeric columns", b"a,b\n1,2\n3,5\n", ["corr", "f.csv"], None, 0),
    ("corr: one numeric column", b"a,b\n1,x\n2,y\n", ["corr", "f.csv"], "stats.corr-columns", 2),
    ("corr: cols names one column", b"a,b\n1,2\n3,5\n", ["corr", "--cols", "a", "f.csv"], "stats.corr-columns", 2),
    ("corr: text column named", b"a,b,c\n1,2,x\n3,5,y\n", ["corr", "--cols", "a,c", "f.csv"], "stats.non-numeric", 2),
    ("corr: empty column named", b"a,b,c\n1,2,\n3,5,\n", ["corr", "--cols", "a,c", "f.csv"], "stats.non-numeric", 2),
    ("corr: no such column", b"a,b\n1,2\n3,5\n", ["corr", "--cols", "a,zz", "f.csv"], "stats.no-such-column", 2),
    ("corr: ambiguous column", b"a,a,b\n1,2,3\n3,5,6\n", ["corr", "--cols", "a,b", "f.csv"], "stats.ambiguous-column", 2),
    ("corr: empty cols item", b"a,b\n1,2\n3,5\n", ["corr", "--cols", "a,,b", "f.csv"], "args.bad-value", 2),
    ("corr: pairs exactly at the limit", b"a,b,c\n1,2,3\n3,5,4\n", ["corr", "--max-pairs", "3", "f.csv"], None, 0),
    ("corr: one pair over the limit", b"a,b,c\n1,2,3\n3,5,4\n", ["corr", "--max-pairs", "2", "f.csv"], "limit.too-many-pairs", 8),
    ("corr: max-pairs zero", b"a,b\n1,2\n3,5\n", ["corr", "--max-pairs", "0", "f.csv"], "args.bad-value", 2),
    ("corr: max-pairs at the ceiling", b"a,b\n1,2\n3,5\n", ["corr", "--max-pairs", "32640", "f.csv"], None, 0),
    ("corr: max-pairs above the ceiling", b"a,b\n1,2\n3,5\n", ["corr", "--max-pairs", "32641", "f.csv"], "args.bad-value", 2),
    ("corr: method unknown", b"a,b\n1,2\n3,5\n", ["corr", "--method", "tau", "f.csv"], "args.bad-value", 2),
    ("corr: kendall", b"a,b\n1,2\n3,5\n", ["corr", "--method", "kendall", "f.csv"], None, 0),
    ("corr: kendall with ties, a missing cell, three columns", b"a,b,c\n1,2,x\n1,,y\n3,2,z\n3,2,w\n", ["corr", "--method", "kendall", "f.csv"], None, 0),
    ("corr: spearman", b"a,b\n1,2\n3,5\n", ["corr", "--method", "spearman", "f.csv"], None, 0),
    ("corr: spearman with ties and a missing cell", b"a,b\n1,2\n1,\n3,2\n", ["corr", "--method", "spearman", "--missing", "listwise", "f.csv"], None, 0),
    ("corr: missing mode unknown", b"a,b\n1,2\n3,5\n", ["corr", "--missing", "mean", "f.csv"], "args.bad-value", 2),
    ("corr: ragged", b"a,b\n1,2\n3\n", ["corr", "f.csv"], "stats.csv-ragged", 8),
    ("corr: no header", b"", ["corr", "f.csv"], "stats.no-header", 8),
    ("corr: rows limit", b"a,b\n1,2\n3,5\n7,8\n", ["corr", "--max-rows", "2", "f.csv"], "limit.too-many-rows", 8),
    ("corr: no file", b"a,b\n1,2\n", ["corr"], "args.missing-operand", 2),
    ("corr: not a verb", b"a,b\n1,2\n", ["correlate", "f.csv"], "stats.bad-verb", 2),
    # test t1 (docs/test.md)
    ("test: t1", b"a\n1\n2\n4\n", ["test", "t1", "--cols", "a", "f.csv"], None, 0),
    ("test: t1 two columns, level and mu0", b"a,b\n1,2\n2,5\n4,1\n", ["test", "t1", "--cols", "a,b", "--mu0", "-1.5", "--level", "0.999999999", "f.csv"], None, 0),
    ("test: t1 with one value, a constant column", b"a,b\n1,3\n,3\n", ["test", "t1", "--cols", "a,b", "f.csv"], None, 0),
    ("test: no kind", b"a\n1\n", ["test"], "args.missing-operand", 2),
    ("test: unknown kind", b"a\n1\n", ["test", "t9", "--cols", "a", "f.csv"], "args.bad-value", 2),
    ("test: no file", b"a\n1\n", ["test", "t1", "--cols", "a"], "args.missing-operand", 2),
    ("test: two files", b"a\n1\n", ["test", "t1", "--cols", "a", "f.csv", "f.csv"], "args.too-many-operands", 2),
    ("test: no columns", b"a\n1\n2\n", ["test", "t1", "f.csv"], "stats.test-columns", 2),
    ("test: text column", b"a,b\n1,x\n2,y\n", ["test", "t1", "--cols", "b", "f.csv"], "stats.non-numeric", 2),
    ("test: empty column", b"a,b\n1,\n2,\n", ["test", "t1", "--cols", "a,b", "f.csv"], "stats.non-numeric", 2),
    ("test: no such column", b"a\n1\n", ["test", "t1", "--cols", "z", "f.csv"], "stats.no-such-column", 2),
    ("test: ambiguous column", b"a,a\n1,2\n", ["test", "t1", "--cols", "a", "f.csv"], "stats.ambiguous-column", 2),
    ("test: empty cols item", b"a\n1\n", ["test", "t1", "--cols", "a,", "f.csv"], "args.bad-value", 2),
    ("test: mu0 not a number", b"a\n1\n", ["test", "t1", "--cols", "a", "--mu0", "abc", "f.csv"], "args.bad-value", 2),
    ("test: mu0 nan", b"a\n1\n", ["test", "t1", "--cols", "a", "--mu0", "nan", "f.csv"], "args.bad-value", 2),
    ("test: mu0 out of range", b"a\n1\n", ["test", "t1", "--cols", "a", "--mu0", "1e999", "f.csv"], "args.bad-value", 2),
    ("test: level zero", b"a\n1\n", ["test", "t1", "--cols", "a", "--level", "0", "f.csv"], "args.bad-value", 2),
    ("test: level one", b"a\n1\n", ["test", "t1", "--cols", "a", "--level", "1", "f.csv"], "args.bad-value", 2),
    ("test: level above one", b"a\n1\n", ["test", "t1", "--cols", "a", "--level", "1.5", "f.csv"], "args.bad-value", 2),
    ("test: level ten decimals", b"a\n1\n", ["test", "t1", "--cols", "a", "--level", "0.1234567891", "f.csv"], "args.bad-value", 2),
    ("test: level text", b"a\n1\n", ["test", "t1", "--cols", "a", "--level", "high", "f.csv"], "args.bad-value", 2),
    ("test: ragged", b"a,b\n1,2\n3\n", ["test", "t1", "--cols", "a", "f.csv"], "stats.csv-ragged", 8),
    ("test: no header", b"", ["test", "t1", "--cols", "a", "f.csv"], "stats.no-header", 8),
    ("test: rows limit", b"a\n1\n2\n3\n", ["test", "t1", "--cols", "a", "--max-rows", "2", "f.csv"], "limit.too-many-rows", 8),
    # test t2 and paired (docs/test.md)
    ("test: t2", b"v,g\n1,a\n2,a\n4,b\n7,b\n", ["test", "t2", "--cols", "v", "--group", "g", "f.csv"], None, 0),
    ("test: t2 pooled, level, two columns", b"v,w,g\n1,5,a\n2,3,a\n4,3,b\n7,1,b\n8,2,b\n", ["test", "t2", "--cols", "v,w", "--group", "g", "--variance", "pooled", "--level", "0.9", "f.csv"], None, 0),
    ("test: t2 groups of one row each", b"v,g\n1,a\n2,b\n", ["test", "t2", "--cols", "v", "--group", "g", "f.csv"], None, 0),
    ("test: t2 constant groups", b"v,g\n3,a\n3,a\n3,b\n3,b\n3,b\n", ["test", "t2", "--cols", "v", "--group", "g", "f.csv"], None, 0),
    ("test: t2 constant groups pooled", b"v,g\n3,a\n3,a\n5,b\n5,b\n", ["test", "t2", "--cols", "v", "--group", "g", "--variance", "pooled", "f.csv"], None, 0),
    ("test: t2 missing group cells and values", b"v,g\n1,a\n,a\n2,\n4,b\n7,b\n9,a\n", ["test", "t2", "--cols", "v", "--group", "g", "f.csv"], None, 0),
    ("test: t2 a group with no value", b"v,g\n1,a\n2,a\n,b\n,b\n", ["test", "t2", "--cols", "v", "--group", "g", "f.csv"], None, 0),
    ("test: t2 no group flag", b"v,g\n1,a\n2,b\n", ["test", "t2", "--cols", "v", "f.csv"], "stats.test-columns", 2),
    ("test: t2 no columns", b"v,g\n1,a\n2,b\n", ["test", "t2", "--group", "g", "f.csv"], "stats.test-columns", 2),
    ("test: t2 one level", b"v,g\n1,a\n2,a\n", ["test", "t2", "--cols", "v", "--group", "g", "f.csv"], "stats.group-levels", 2),
    ("test: t2 three levels", b"v,g\n1,a\n2,b\n3,c\n", ["test", "t2", "--cols", "v", "--group", "g", "f.csv"], "stats.group-levels", 2),
    ("test: t2 group column empty", b"v,g\n1,\n2,\n", ["test", "t2", "--cols", "v", "--group", "g", "f.csv"], "stats.group-levels", 2),
    ("test: t2 no such group column", b"v,g\n1,a\n2,b\n", ["test", "t2", "--cols", "v", "--group", "z", "f.csv"], "stats.no-such-column", 2),
    ("test: t2 ambiguous group column", b"v,g,g\n1,a,a\n2,b,b\n", ["test", "t2", "--cols", "v", "--group", "g", "f.csv"], "stats.ambiguous-column", 2),
    ("test: t2 no such column", b"v,g\n1,a\n2,b\n", ["test", "t2", "--cols", "z", "--group", "g", "f.csv"], "stats.no-such-column", 2),
    ("test: t2 text column", b"v,w,g\n1,x,a\n2,y,b\n", ["test", "t2", "--cols", "w", "--group", "g", "f.csv"], "stats.non-numeric", 2),
    ("test: t2 text in one group", b"v,g\n1,a\n2,a\nx,b\n3,b\n", ["test", "t2", "--cols", "v", "--group", "g", "f.csv"], "stats.non-numeric", 2),
    ("test: t2 variance unknown", b"v,g\n1,a\n2,b\n", ["test", "t2", "--cols", "v", "--group", "g", "--variance", "equal", "f.csv"], "args.bad-value", 2),
    ("test: t2 groups limit", b"v,g\n1,a\n2,b\n", ["test", "t2", "--cols", "v", "--group", "g", "--max-groups", "1", "f.csv"], "limit.too-many-groups", 8),
    ("test: t2 rows limit", b"v,g\n1,a\n2,b\n3,b\n", ["test", "t2", "--cols", "v", "--group", "g", "--max-rows", "2", "f.csv"], "limit.too-many-rows", 8),
    ("test: t2 ragged", b"v,g\n1,a\n2\n", ["test", "t2", "--cols", "v", "--group", "g", "f.csv"], "stats.csv-ragged", 8),
    ("test: t2 level out of range", b"v,g\n1,a\n2,b\n", ["test", "t2", "--cols", "v", "--group", "g", "--level", "1", "f.csv"], "args.bad-value", 2),
    ("test: paired", b"a,b\n1,2\n2,5\n4,1\n", ["test", "paired", "--pairs", "a:b", "f.csv"], None, 0),
    ("test: paired two pairs, one reversed, level", b"a,b,c\n1,2,7\n2,5,3\n4,1,2\n", ["test", "paired", "--pairs", "a:b,c:a,b:c", "--level", "0.99", "f.csv"], None, 0),
    ("test: paired with missing cells", b"a,b\n1,\n,2\n3,4\n4,7\n5,5.5\n", ["test", "paired", "--pairs", "a:b", "f.csv"], None, 0),
    ("test: paired a column with itself", b"a,b\n1,2\n2,5\n4,1\n", ["test", "paired", "--pairs", "a:a", "f.csv"], None, 0),
    ("test: paired one pair of rows", b"a,b\n1,2\n", ["test", "paired", "--pairs", "a:b", "f.csv"], None, 0),
    ("test: paired no pairs", b"a,b\n1,2\n", ["test", "paired", "f.csv"], "stats.test-columns", 2),
    ("test: paired --cols is not enough", b"a,b\n1,2\n", ["test", "paired", "--cols", "a,b", "f.csv"], "stats.test-columns", 2),
    ("test: paired item without colon", b"a,b\n1,2\n", ["test", "paired", "--pairs", "a", "f.csv"], "args.bad-value", 2),
    ("test: paired empty name", b"a,b\n1,2\n", ["test", "paired", "--pairs", "a:", "f.csv"], "args.bad-value", 2),
    ("test: paired empty first name", b"a,b\n1,2\n", ["test", "paired", "--pairs", ":b", "f.csv"], "args.bad-value", 2),
    ("test: paired empty item", b"a,b\n1,2\n", ["test", "paired", "--pairs", "a:b,", "f.csv"], "args.bad-value", 2),
    ("test: paired no such column", b"a,b\n1,2\n", ["test", "paired", "--pairs", "a:z", "f.csv"], "stats.no-such-column", 2),
    ("test: paired ambiguous column", b"a,a\n1,2\n", ["test", "paired", "--pairs", "a:a", "f.csv"], "stats.ambiguous-column", 2),
    ("test: paired second name ambiguous", b"a,b,b\n1,2,3\n", ["test", "paired", "--pairs", "a:b", "f.csv"], "stats.ambiguous-column", 2),
    ("test: paired second name missing", b"a,b\n1,2\n", ["test", "paired", "--pairs", "b:z", "f.csv"], "stats.no-such-column", 2),
    ("test: paired pairs at the limit", b"a,b,c\n1,2,3\n2,5,1\n", ["test", "paired", "--pairs", "a:b", "--max-pairs", "1", "f.csv"], None, 0),
    # test chisq (docs/test.md section 10)
    ("test: chisq independence", b"a,g\nx,u\nx,v\ny,u\ny,v\ny,v\n", ["test", "chisq", "--cols", "a", "--vs", "g", "f.csv"], None, 0),
    ("test: chisq two columns, numbers as levels", b"a,b,g\n1,x,u\n10,y,v\n9,x,u\n1,y,v\n", ["test", "chisq", "--cols", "a,b", "--vs", "g", "f.csv"], None, 0),
    ("test: chisq fit", b"a\nx\nx\ny\nz\n", ["test", "chisq", "--cols", "a", "--expected", "0.5,0.25,0.25", "f.csv"], None, 0),
    ("test: chisq fit, nine decimals", b"a\nx\ny\ny\n", ["test", "chisq", "--cols", "a", "--expected", "0.333333333,0.666666667", "f.csv"], None, 0),
    ("test: chisq missing cells", b"a,g\nx,u\n,v\ny,\ny,v\nx,v\n", ["test", "chisq", "--cols", "a", "--vs", "g", "f.csv"], None, 0),
    ("test: chisq one level", b"a,g\nx,u\nx,v\n", ["test", "chisq", "--cols", "a", "--vs", "g", "f.csv"], None, 0),
    ("test: chisq no complete row", b"a,g\nx,\n,v\n", ["test", "chisq", "--cols", "a", "--vs", "g", "f.csv"], None, 0),
    ("test: chisq fit with one level", b"a\nx\nx\n", ["test", "chisq", "--cols", "a", "--expected", "1", "f.csv"], None, 0),
    ("test: chisq neither vs nor expected", b"a,g\nx,u\n", ["test", "chisq", "--cols", "a", "f.csv"], "stats.test-columns", 2),
    ("test: chisq both vs and expected", b"a,g\nx,u\n", ["test", "chisq", "--cols", "a", "--vs", "g", "--expected", "1", "f.csv"], "stats.test-columns", 2),
    ("test: chisq no columns", b"a,g\nx,u\n", ["test", "chisq", "--vs", "g", "f.csv"], "stats.test-columns", 2),
    ("test: chisq empty cols item", b"a,g\nx,u\n", ["test", "chisq", "--cols", "a,", "--vs", "g", "f.csv"], "args.bad-value", 2),
    ("test: chisq no such column", b"a,g\nx,u\n", ["test", "chisq", "--cols", "z", "--vs", "g", "f.csv"], "stats.no-such-column", 2),
    ("test: chisq no such vs column", b"a,g\nx,u\n", ["test", "chisq", "--cols", "a", "--vs", "z", "f.csv"], "stats.no-such-column", 2),
    ("test: chisq ambiguous column", b"a,a,g\nx,x,u\n", ["test", "chisq", "--cols", "a", "--vs", "g", "f.csv"], "stats.ambiguous-column", 2),
    ("test: chisq ambiguous vs column", b"a,g,g\nx,u,u\n", ["test", "chisq", "--cols", "a", "--vs", "g", "f.csv"], "stats.ambiguous-column", 2),
    ("test: chisq expected not decimals", b"a\nx\ny\n", ["test", "chisq", "--cols", "a", "--expected", "a,b", "f.csv"], "stats.expected-proportions", 2),
    ("test: chisq expected empty item", b"a\nx\ny\n", ["test", "chisq", "--cols", "a", "--expected", "0.5,0.5,", "f.csv"], "stats.expected-proportions", 2),
    ("test: chisq expected sum", b"a\nx\ny\n", ["test", "chisq", "--cols", "a", "--expected", "0.5,0.6", "f.csv"], "stats.expected-proportions", 2),
    ("test: chisq expected sum below one", b"a\nx\ny\n", ["test", "chisq", "--cols", "a", "--expected", "0.4,0.4", "f.csv"], "stats.expected-proportions", 2),
    ("test: chisq expected two points", b"a\nx\ny\n", ["test", "chisq", "--cols", "a", "--expected", "0.5.5,0.45", "f.csv"], "stats.expected-proportions", 2),
    ("test: chisq expected zero", b"a\nx\ny\n", ["test", "chisq", "--cols", "a", "--expected", "0,1", "f.csv"], "stats.expected-proportions", 2),
    ("test: chisq expected ten decimals", b"a\nx\ny\n", ["test", "chisq", "--cols", "a", "--expected", "0.0000000001,0.999999999", "f.csv"], "stats.expected-proportions", 2),
    ("test: chisq expected negative", b"a\nx\ny\n", ["test", "chisq", "--cols", "a", "--expected", "1.5,-0.5", "f.csv"], "stats.expected-proportions", 2),
    ("test: chisq expected count differs", b"a\nx\ny\nz\n", ["test", "chisq", "--cols", "a", "--expected", "0.5,0.5", "f.csv"], "stats.expected-proportions", 2),
    ("test: chisq expected count above the levels", b"a\nx\ny\n", ["test", "chisq", "--cols", "a", "--expected", "0.2,0.3,0.5", "f.csv"], "stats.expected-proportions", 2),
    ("test: chisq groups limit", b"a,g\nx,u\ny,v\n", ["test", "chisq", "--cols", "a", "--vs", "g", "--max-groups", "1", "f.csv"], "limit.too-many-groups", 8),
    ("test: chisq groups limit on the second column", b"a,g\nx,u\nx,v\n", ["test", "chisq", "--cols", "a", "--vs", "g", "--max-groups", "1", "f.csv"], "limit.too-many-groups", 8),
    ("test: chisq groups at the limit", b"a,g\nx,u\nx,v\n", ["test", "chisq", "--cols", "a", "--vs", "g", "--max-groups", "2", "f.csv"], None, 0),
    ("test: chisq rows limit", b"a,g\nx,u\ny,v\nx,v\n", ["test", "chisq", "--cols", "a", "--vs", "g", "--max-rows", "2", "f.csv"], "limit.too-many-rows", 8),
    ("test: chisq ragged", b"a,g\nx,u\ny\n", ["test", "chisq", "--cols", "a", "--vs", "g", "f.csv"], "stats.csv-ragged", 8),
    ("test: chisq no header", b"", ["test", "chisq", "--cols", "a", "--vs", "g", "f.csv"], "stats.no-header", 8),
    # test prop (docs/test.md section 11)
    ("test: prop one sample", b"o\ny\nn\ny\ny\n", ["test", "prop", "--cols", "o", "--success", "y", "f.csv"], None, 0),
    ("test: prop p0, level, two columns", b"o,w\ny,n\nn,n\ny,y\n,y\n", ["test", "prop", "--cols", "o,w", "--success", "y", "--p0", "0.333333333", "--level", "0.9", "f.csv"], None, 0),
    ("test: prop two samples", b"o,g\ny,a\ny,a\nn,a\nn,b\ny,b\nn,b\n", ["test", "prop", "--cols", "o", "--success", "y", "--group", "g", "f.csv"], None, 0),
    ("test: prop no success in the file", b"o\nn\nn\n", ["test", "prop", "--cols", "o", "--success", "y", "f.csv"], None, 0),
    ("test: prop every row a success", b"o,g\ny,a\ny,b\n", ["test", "prop", "--cols", "o", "--success", "y", "--group", "g", "f.csv"], None, 0),
    ("test: prop no complete row", b"o,g\n,a\ny,\n", ["test", "prop", "--cols", "o", "--success", "y", "f.csv"], None, 0),
    ("test: prop a group without a trial", b"o,g\ny,a\n,b\nn,a\n", ["test", "prop", "--cols", "o", "--success", "y", "--group", "g", "f.csv"], None, 0),
    ("test: prop no success flag", b"o\ny\n", ["test", "prop", "--cols", "o", "f.csv"], "stats.test-columns", 2),
    ("test: prop no columns", b"o\ny\n", ["test", "prop", "--success", "y", "f.csv"], "stats.test-columns", 2),
    ("test: prop empty cols item", b"o\ny\n", ["test", "prop", "--cols", "o,", "--success", "y", "f.csv"], "args.bad-value", 2),
    ("test: prop p0 zero", b"o\ny\n", ["test", "prop", "--cols", "o", "--success", "y", "--p0", "0", "f.csv"], "args.bad-value", 2),
    ("test: prop p0 one", b"o\ny\n", ["test", "prop", "--cols", "o", "--success", "y", "--p0", "1", "f.csv"], "args.bad-value", 2),
    ("test: prop p0 above one", b"o\ny\n", ["test", "prop", "--cols", "o", "--success", "y", "--p0", "1.5", "f.csv"], "args.bad-value", 2),
    ("test: prop p0 text", b"o\ny\n", ["test", "prop", "--cols", "o", "--success", "y", "--p0", "half", "f.csv"], "args.bad-value", 2),
    ("test: prop p0 ten decimals", b"o\ny\n", ["test", "prop", "--cols", "o", "--success", "y", "--p0", "0.1234567891", "f.csv"], "args.bad-value", 2),
    ("test: prop p0 with group", b"o,g\ny,a\nn,b\n", ["test", "prop", "--cols", "o", "--success", "y", "--group", "g", "--p0", "0.5", "f.csv"], "args.bad-value", 2),
    ("test: prop one group level", b"o,g\ny,a\nn,a\n", ["test", "prop", "--cols", "o", "--success", "y", "--group", "g", "f.csv"], "stats.group-levels", 2),
    ("test: prop three group levels", b"o,g\ny,a\nn,b\ny,c\n", ["test", "prop", "--cols", "o", "--success", "y", "--group", "g", "f.csv"], "stats.group-levels", 2),
    ("test: prop no such column", b"o\ny\n", ["test", "prop", "--cols", "z", "--success", "y", "f.csv"], "stats.no-such-column", 2),
    ("test: prop no such group column", b"o\ny\n", ["test", "prop", "--cols", "o", "--success", "y", "--group", "z", "f.csv"], "stats.no-such-column", 2),
    ("test: prop ambiguous column", b"o,o\ny,n\n", ["test", "prop", "--cols", "o", "--success", "y", "f.csv"], "stats.ambiguous-column", 2),
    ("test: prop levels limit", b"o\ny\nn\n", ["test", "prop", "--cols", "o", "--success", "y", "--max-groups", "1", "f.csv"], "limit.too-many-groups", 8),
    ("test: prop rows limit", b"o\ny\nn\ny\n", ["test", "prop", "--cols", "o", "--success", "y", "--max-rows", "2", "f.csv"], "limit.too-many-rows", 8),
    ("test: prop ragged", b"o,g\ny,a\nn\n", ["test", "prop", "--cols", "o", "--success", "y", "--group", "g", "f.csv"], "stats.csv-ragged", 8),
    ("test: paired text column", b"a,b\n1,x\n2,y\n", ["test", "paired", "--pairs", "a:b", "f.csv"], "stats.non-numeric", 2),
    ("test: paired pairs limit", b"a,b,c\n1,2,3\n", ["test", "paired", "--pairs", "a:b,b:c", "--max-pairs", "1", "f.csv"], "limit.too-many-pairs", 8),
    ("test: paired rows limit", b"a,b\n1,2\n2,3\n3,4\n", ["test", "paired", "--pairs", "a:b", "--max-rows", "2", "f.csv"], "limit.too-many-rows", 8),
    ("test: paired ragged", b"a,b\n1,2\n3\n", ["test", "paired", "--pairs", "a:b", "f.csv"], "stats.csv-ragged", 8),
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
            k = int(r.next() % 30)
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
            elif k == 10:
                argv = ["corr", "f.csv"]
            elif k == 11:
                argv = ["corr", "--missing", "listwise", "f.csv"]
            elif k == 12:
                argv = ["corr", "--max-pairs", "1", "--max-rows", "3", "f.csv"]
            elif k == 13:
                argv = ["corr", "--cols", "a,b", "--format", "text", "f.csv"]
            elif k == 14:
                argv = ["corr", "--method", "spearman", "f.csv"]
            elif k == 15:
                argv = ["corr", "--method", "spearman", "--missing", "listwise", "f.csv"]
            elif k == 16:
                argv = ["corr", "--method", "kendall", "f.csv"]
            elif k == 17:
                argv = ["corr", "--method", "kendall", "--missing", "listwise", "--max-pairs", "1", "f.csv"]
            elif k == 18:
                argv = ["test", "t1", "--cols", "a", "f.csv"]
            elif k == 19:
                argv = ["test", "t1", "--cols", "a,b", "--mu0", "1", "--level", "0.9", "f.csv"]
            elif k == 20:
                argv = ["test", "t2", "--cols", "a", "--group", "b", "f.csv"]
            elif k == 21:
                argv = ["test", "t2", "--cols", "a,b", "--group", "a", "--variance", "pooled", "--max-groups", "3", "f.csv"]
            elif k == 22:
                argv = ["test", "paired", "--pairs", "a:b", "f.csv"]
            elif k == 27:
                argv = ["test", "prop", "--cols", "a", "--success", "1", "f.csv"]
            elif k == 28:
                argv = ["test", "prop", "--cols", "a,b", "--success", "a", "--group", "b", "--level", "0.99", "f.csv"]
            elif k == 29:
                argv = ["test", "prop", "--cols", "b", "--success", "", "--p0", "0.000000001", "--max-rows", "5", "f.csv"]
            elif k == 24:
                argv = ["test", "chisq", "--cols", "a", "--vs", "b", "f.csv"]
            elif k == 25:
                argv = ["test", "chisq", "--cols", "a,b", "--expected", "0.25,0.75", "--max-groups", "3", "f.csv"]
            elif k == 26:
                argv = ["test", "chisq", "--cols", "b,a", "--vs", "a", "--max-rows", "6", "f.csv"]
            elif k == 23:
                argv = ["test", "paired", "--pairs", "b:a,a:a", "--level", "0.9", "--max-rows", "6", "f.csv"]
            try:
                p = call(a.bin, tmp, argv)
            except subprocess.TimeoutExpired:
                bad.append("fuzz %d: hang on %r" % (i, data[:60]))
                continue
            if k == 9 or k == 13:
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
