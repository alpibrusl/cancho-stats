#!/usr/bin/env python3
"""Fills the runs of README.md: every `<!--run:ID-->` ... `<!--/run-->` region is rewritten with the real output of build/stats on the files of docs/, so no output in the README is typed by hand.

    python3 scripts/readme_runs.py [--bin build/stats] [--check]

`--check` fails when the README differs from what the binary prints now (CI runs it on every platform that builds the binary)."""
import argparse
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
BAKERY = "docs/stories/bakery.csv"
READINGS = "docs/example/readings.csv"


def run(binary, args):
    p = subprocess.run([binary] + args, cwd=ROOT, capture_output=True, text=True, timeout=300)
    return p.stdout, p.returncode


def shown(args):
    return "stats " + " ".join(a if re.fullmatch(r"[\w.,:=/-]+", a) else "'%s'" % a for a in args)


def text_run(binary, args, lines):
    out, rc = run(binary, args)
    if rc != 0:
        raise SystemExit("readme_runs: `%s` exited %d" % (shown(args), rc))
    rows = out.rstrip("\n").split("\n")
    body = "\n".join(rows[:lines]) + ("\n..." if len(rows) > lines else "")
    return "$ %s\n%s" % (shown(args), body)


def refusal_run(binary, args):
    out, rc = run(binary, args)
    err = json.loads(out)["error"]
    return "$ %s\n%s\n$ echo $?\n%d" % (shown(args), json.dumps({k: err[k] for k in ("code", "rule", "hint")}, indent=2), rc)


def regions(binary):
    text = ["--format", "text"]
    flow = [["describe", "--cols", "drift"] + text + [READINGS]]
    seedless = ["resample", "perm", "--cols", "bake_min", "--group", "oven", BAKERY]
    seeded = ["resample", "perm", "--cols", "bake_min", "--group", "oven", "--seed", "1", "--permutations", "9999"] + text + [BAKERY]
    return {
        "flow": text_run(binary, flow[0], 12) + "\n\n" + refusal_run(binary, seedless),
        "t2": text_run(binary, ["test", "t2", "--cols", "bake_min", "--group", "oven"] + text + [BAKERY], 12),
        "ci": text_run(binary, ["ci", "diff", "--cols", "bake_min", "--group", "oven"] + text + [BAKERY], 10),
        "perm": text_run(binary, seeded, 12),
        "agent": ("# 1. the agent asks the tool what it is, not this file\n" + text_run(binary, ["skill"], 4) + "\n\n# 2. it runs the test it wants, without a seed: a refusal, as data\n" + refusal_run(binary, seedless)
                  + "\n\n# 3. it applies the hint (a seed) and runs again\n" + text_run(binary, seeded, 9)),
    }


def render(readme, binary):
    runs = regions(binary)

    def sub(m):
        if m.group(1) not in runs:
            raise SystemExit("readme_runs: no run named %r" % m.group(1))
        return "<!--run:%s-->\n```console\n%s\n```\n<!--/run-->" % (m.group(1), runs[m.group(1)])
    return re.sub(r"<!--run:(\w+)-->.*?<!--/run-->", sub, readme, flags=re.S)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    path = ROOT / "README.md"
    old = path.read_text()
    new = render(old, a.bin)
    if a.check:
        if old != new:
            print("readme_runs: README.md is stale; run python3 scripts/readme_runs.py")
            return 1
        print("readme_runs: up to date")
        return 0
    path.write_text(new)
    print("readme_runs: wrote README.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
