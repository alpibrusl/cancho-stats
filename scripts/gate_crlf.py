#!/usr/bin/env python3
"""A CRLF file is the same table as its LF twin (RFC 4180: CRLF ends a record, so the carriage return is not part of the
last cell). Regression gate for cancho-stats#62: `describe`, `corr` and `test` must give the same bytes on both.

    python3 scripts/gate_crlf.py [--bin build/stats]

Exit 0: all of it holds. Exit 1: a violation. Exit 2: the harness failed.
"""
import argparse
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
LF = b"g,x,y\nA,1,2.5\nA,2,3.5\nA,4,1\nB,3,9\nB,5,8\nB,6,7.5\n"
CASES = [
    ["describe", "--cols", "x,y", "f.csv"],
    ["describe", "--cols", "y", "--by", "g", "f.csv"],
    ["corr", "--cols", "x,y", "f.csv"],
    ["test", "t2", "--cols", "y", "--group", "g", "f.csv"],
    ["test", "mannwhitney", "--cols", "y", "--group", "g", "f.csv"],
]


def run(binary, args, data):
    with tempfile.TemporaryDirectory() as d:
        (pathlib.Path(d) / "f.csv").write_bytes(data)
        p = subprocess.run([binary, "--root", d] + args, capture_output=True, timeout=60)
        return p.returncode, p.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    a = ap.parse_args()
    crlf = LF.replace(b"\n", b"\r\n")
    bad = 0
    for args in CASES:
        r1, o1 = run(a.bin, args, LF)
        r2, o2 = run(a.bin, args, crlf)
        if r1 != 0 or r2 != 0 or o1 != o2:
            bad += 1
            print("FAIL", " ".join(args), "lf", r1, "crlf", r2, "same bytes" if o1 == o2 else "different bytes")
    print("crlf: %d cases, %d violations" % (len(CASES), bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
