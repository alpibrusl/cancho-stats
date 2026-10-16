#!/usr/bin/env python3
"""Rewrite the replay of docs/index.html from REAL runs of build/stats (the memory of this project: the page shows the real CLI and is updated when the CLI grows).

    python3 scripts/update_replay.py            # rewrite the `var lines = [...]` block of docs/index.html
    python3 scripts/update_replay.py --check    # exit 1 if the page differs from what the binary prints now

The replay is `cat readings.csv`, `stats describe readings.csv` on docs/example/readings.csv, and `stats describe broken.csv` on a file with a short record (its exit code is printed).
"""
import json
import pathlib
import re
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
BROKEN = "city,temp\nA1,21.4\nA2\nA3,22.0\n"


def lines():
    csv = (ROOT / "docs" / "example" / "readings.csv").read_text()
    out = [["d", "$ cat readings.csv"]] + [["a", ln] for ln in csv.rstrip("\n").split("\n")]
    binary = str(ROOT / "build" / "stats")
    ok = subprocess.run([binary, "describe", "--root", str(ROOT / "docs" / "example"), "readings.csv"], capture_output=True, text=True, check=True)
    out += [["d", "$ stats describe readings.csv"], ["g", ok.stdout.rstrip("\n")]]
    with tempfile.TemporaryDirectory() as t:
        (pathlib.Path(t) / "broken.csv").write_text(BROKEN)
        bad = subprocess.run([binary, "describe", "--root", t, "broken.csv"], capture_output=True, text=True)
    out += [["d", '$ stats describe broken.csv; echo "exit $?"'], ["b", bad.stdout.rstrip("\n")], ["b", "exit %d" % bad.returncode]]
    return out


def render(items):
    return "var lines = " + json.dumps(items, indent=2) + ";"


def main():
    page = ROOT / "docs" / "index.html"
    text = page.read_text()
    m = re.search(r"var lines = \[.*?\n\];", text, re.S)
    if not m:
        print("update_replay: no `var lines = [...]` block in docs/index.html")
        return 2
    new = text[: m.start()] + render(lines()) + text[m.end():]
    if "--check" in sys.argv:
        if new != text:
            print("update_replay: docs/index.html shows other bytes than build/stats prints now")
            return 1
        print("update_replay: the replay equals the real runs")
        return 0
    page.write_text(new)
    print("update_replay: wrote the replay (%d lines)" % len(lines()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
