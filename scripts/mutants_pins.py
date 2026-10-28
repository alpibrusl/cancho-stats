#!/usr/bin/env python3
"""Mutants of scripts/check_pins.sh (docs/numerics.md section 2.1): each is applied to a COPY of the files the check reads and must make it exit 1 with the
message of the id it is for. The unmutated copy is the control and must exit 0, or every kill below would be a broken harness.

    python3 scripts/mutants_pins.py            # run; exit 1 if a mutant survives or is killed for another reason than its id
    python3 scripts/mutants_pins.py --check    # only verify that every pattern applies exactly once
"""
import os, pathlib, re, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
TOML = "cancho.toml"
DOC = "docs/numerics.md"
FIRST = "[dependencies.ml-rng]\ngit = \"https://github.com/alpibrusl/cancho-ml\"\nrev = \"99e9b64afb4f1ca191a4209dea1306720da43b2f\"\npath = \".cancho-vcs/ml.rng\"\n"
REV_ML = "99e9b64afb4f1ca191a4209dea1306720da43b2f"
REV_TAB = "202b214372a03ec6ed8ca841357f3d2ed4de8aa0"
FLT = "[dependencies.flt]\ngit = \"https://github.com/alpibrusl/cancho-table\"\nrev = \"%s\"\npath = \".cancho-vcs/flt\"\n" % REV_TAB

# (name, file, old, new, id) ; special mutants (no text edit) have file None and a callable
MUTANTS = [
    ("block-removed", TOML, FIRST, "", "PN1"),
    ("unlisted-block", TOML, FLT, FLT + "\n[dependencies.query]\ngit = \"https://github.com/alpibrusl/cancho-table\"\nrev = \"%s\"\npath = \".cancho-vcs/query\"\n" % REV_TAB, "PN1"),
    ("block-duplicated", TOML, FLT, FLT + "\n" + FLT, "PN2"),
    ("same-store-other-label", TOML, FLT, FLT + "\n[dependencies.flt2]\ngit = \"https://github.com/alpibrusl/cancho-table\"\nrev = \"%s\"\npath = \".cancho-vcs/flt\"\n" % REV_TAB, "PN2"),
    ("rev-39-digits", TOML, "rev = \"%s\"\npath = \".cancho-vcs/flt\"" % REV_TAB, "rev = \"%s\"\npath = \".cancho-vcs/flt\"" % REV_TAB[:39], "PN3"),
    ("rev-not-hex", TOML, "rev = \"%s\"\npath = \".cancho-vcs/flt\"" % REV_TAB, "rev = \"%sz\"\npath = \".cancho-vcs/flt\"" % REV_TAB[:39], "PN3"),
    ("rev-differs-in-one-block", TOML, "rev = \"%s\"\npath = \".cancho-vcs/flt\"" % REV_TAB, "rev = \"%s1\"\npath = \".cancho-vcs/flt\"" % REV_TAB[:39], "PN3"),
    ("rev-is-not-the-documented", DOC, "cancho-table rev: `%s`" % REV_TAB, "cancho-table rev: `%s`" % ("0" * 40), "PN3"),
    ("document-has-no-rev-line", DOC, "cancho-ml rev: `%s`" % REV_ML, "cancho-ml rev is somewhere else", "PN3"),
    ("compiler-rev-moved", TOML, "cancho = \"f4ead824f33b291a0ad7620a34a0c1346a7abe4f\"", "cancho = \"e47d647e47d647e47d647e47d647e47d647e47d6\"", "PN3"),
    ("git-elsewhere", TOML, "git = \"https://github.com/alpibrusl/cancho-table\"\nrev = \"%s\"\npath = \".cancho-vcs/flt\"" % REV_TAB, "git = \"https://github.com/someone/cancho-table\"\nrev = \"%s\"\npath = \".cancho-vcs/flt\"" % REV_TAB, "PN4"),
    ("ml-git-elsewhere", TOML, FIRST, FIRST.replace("alpibrusl/cancho-ml", "someone/cancho-ml"), "PN4"),
    ("path-of-another-module", TOML, "path = \".cancho-vcs/flt\"", "path = \".cancho-vcs/pow\"", "PN5"),
    ("path-without-prefix", TOML, "path = \".cancho-vcs/flt\"", "path = \"flt\"", "PN5"),
    ("sources-leave-the-tree", TOML, "[package]\n", "[[bin]]\nname = \"x\"\nsources = [\"../cancho-table/tools/table\"]\n\n[package]\n", "PN6"),
    ("set-differs-in-document", DOC, "| cancho-table | flt | flt |\n", "", "PN7"),
    ("document-lists-an-extra", DOC, "| cancho-table | flt | flt |\n", "| cancho-table | flt | flt |\n| cancho-table | query | query |\n", "PN7"),
    ("document-has-no-set", DOC, "<!-- pins:begin -->", "<!-- pins:begun -->", "PN7"),
]
SPECIAL = ["copy-in-lib", "copy-in-tools", "deps-missing-module", "deps-complete-control"]


def tree(tmp):
    for f in (TOML, DOC, "scripts/check_pins.sh"):
        d = tmp / f
        d.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / f, d)
    (tmp / "lib").mkdir(exist_ok=True)
    (tmp / "tools").mkdir(exist_ok=True)


def run(tmp):
    r = subprocess.run(["sh", "scripts/check_pins.sh"], cwd=tmp, capture_output=True, text=True)
    return r.returncode, r.stderr + r.stdout


def fake_deps(tmp, drop=None):
    mods = re.findall(r"^\| \S+ \| \S+ \| (\S+) \|$", (tmp / DOC).read_text(), re.M)
    d = tmp / "build" / "deps"
    d.mkdir(parents=True)
    for i, m in enumerate(mods):
        if m != drop:
            (d / ("%02d.cho" % i)).write_text("module %s;\n" % m)


def problems():
    out = []
    names = set()
    for name, file, old, new, pid in MUTANTS:
        if name in names:
            out.append("mutant %s: named twice" % name)
        names.add(name)
        t = (ROOT / file).read_text()
        n = t.count(old)
        if n != 1:
            out.append("mutant %s: pattern occurs %d times in %s (must be exactly 1)" % (name, n, file))
        if old == new:
            out.append("mutant %s: replacement equals pattern" % name)
    return out


def main():
    bad = problems()
    if bad:
        print("\n".join(bad))
        return 1
    if "--check" in sys.argv:
        print("mutants_pins: %d text mutants, every pattern applies once" % len(MUTANTS))
        return 0
    survivors, killed = [], 0
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        tree(tmp)
        code, msg = run(tmp)
        if code != 0:
            print("mutants_pins: the unmutated control FAILS: " + msg)
            return 1
        print("control (unmutated) passes")
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        tree(tmp)
        fake_deps(tmp)
        code, msg = run(tmp)
        if code != 0:
            print("mutants_pins: the control with a complete build/deps FAILS: " + msg)
            return 1
        print("control (complete build/deps) passes")
    cases = [(n, f, o, w, p) for (n, f, o, w, p) in MUTANTS]
    for name, file, old, new, pid in cases:
        with tempfile.TemporaryDirectory() as t:
            tmp = pathlib.Path(t)
            tree(tmp)
            f = tmp / file
            f.write_text(f.read_text().replace(old, new, 1))
            code, msg = run(tmp)
        ok = code == 1 and ("check_pins: %s:" % pid) in msg
        if ok:
            killed += 1
            print("%-32s killed by %s" % (name, pid))
        else:
            survivors.append(name)
            print("%-32s NOT KILLED BY %s (exit %d): %s" % (name, pid, code, msg.strip()[:160]))
    # special mutants: a copy of a pinned module in the tree; build/deps lacking a module
    for name, pid, where in (("copy-in-lib", "PN6", "lib/stats_copy.cho"), ("copy-in-tools", "PN6", "tools/x/copy.cho")):
        with tempfile.TemporaryDirectory() as t:
            tmp = pathlib.Path(t)
            tree(tmp)
            p = tmp / where
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("module facc;\n\npub fn width() -> [] int {\n    return 73;\n}\n")
            code, msg = run(tmp)
        ok = code == 1 and ("check_pins: %s:" % pid) in msg
        killed += ok
        if not ok:
            survivors.append(name)
        print("%-32s %s" % (name, "killed by " + pid if ok else "NOT KILLED: " + msg.strip()[:160]))
    with tempfile.TemporaryDirectory() as t:
        tmp = pathlib.Path(t)
        tree(tmp)
        fake_deps(tmp, drop="linalg.chol")
        code, msg = run(tmp)
    ok = code == 1 and "check_pins: PN8:" in msg
    killed += ok
    if not ok:
        survivors.append("deps-missing-module")
    print("%-32s %s" % ("deps-missing-module", "killed by PN8" if ok else "NOT KILLED: " + msg.strip()[:160]))
    total = len(cases) + 3
    print("mutants_pins: %d killed of %d, %d surviving" % (killed, total, len(survivors)))
    if survivors:
        print("SURVIVORS: " + ", ".join(survivors))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
