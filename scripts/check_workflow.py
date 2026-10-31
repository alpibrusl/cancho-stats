#!/usr/bin/env python3
"""Structural checks of .github/workflows/ci.yml (docs/numerics.md section 7.4, gate G7), for a machine without `actionlint`. Nothing here runs the workflow.

    python3 scripts/check_workflow.py [FILE]       # exit 1 with the reasons when one of the rules below is broken
    python3 scripts/check_workflow.py --mutants    # each rule is shown able to fail: a copy of the workflow with one fault must be refused

  W1  the top-level `permissions` is exactly {contents: read}
  W2  no job has a `permissions` that grants write (a job's own block, if any, grants read or none only)
  W3  every `uses:` is an action pinned to a 40-hex commit SHA (a local `./` action is allowed), with a version comment on the line
  W4  no `pull_request_target`, and no `secrets.` at all (the workflow needs none)
  W5  the compiler revision is read from cancho.toml with sed, and no 40-hex string is a literal in the workflow other than an action pin
  W6  every script the workflow runs (`scripts/...`, `ci/...`) is a file of the repository
  W7  triggers: push, pull_request and workflow_dispatch; `full-gate` runs only for workflow_dispatch and schedule
  W8  the matrix of `build-test` has a Linux x86-64, a Linux aarch64 and a macOS runner
  W9  every job has `timeout-minutes`
  W10 every runner image is named: no `-latest` anywhere, and every `runs-on` is a name with a version or a matrix expression
"""
import pathlib
import re
import sys

try:
    import yaml
except ImportError:
    print("check_workflow: PyYAML is not installed (pip install -r ci/requirements.txt)", file=sys.stderr)
    sys.exit(2)

ROOT = pathlib.Path(__file__).resolve().parent.parent
DEFAULT = ROOT / ".github" / "workflows" / "ci.yml"
SHA = re.compile(r"^[^@\s]+@[0-9a-f]{40}$")


def load(text):
    d = yaml.safe_load(text)
    if not isinstance(d, dict):
        raise ValueError("not a mapping")
    if True in d and "on" not in d:  # YAML 1.1 reads the key `on` as true
        d["on"] = d.pop(True)
    return d


def check(text, root=ROOT):
    """-> list of (rule, reason)"""
    out = []
    try:
        d = load(text)
    except Exception as e:
        return [("W0", "the workflow does not parse: %s" % e)]
    jobs = d.get("jobs") or {}
    if d.get("permissions") != {"contents": "read"}:
        out.append(("W1", "top-level permissions is %r, not {contents: read}" % (d.get("permissions"),)))
    for name, job in jobs.items():
        p = job.get("permissions")
        if p is not None and (not isinstance(p, dict) or any(v not in ("read", "none") for v in p.values())):
            out.append(("W2", "job %s grants %r; nothing here may write" % (name, p)))
        if "timeout-minutes" not in job:
            out.append(("W9", "job %s has no timeout-minutes" % name))
        runs = str(job.get("runs-on", ""))
        if not runs or not (re.match(r"^[a-z]+-[0-9]+(\.[0-9]+)?(-arm)?$", runs) or "${{ matrix.os }}" in runs):
            out.append(("W10", "job %s runs-on %r is not a named image" % (name, runs)))
        for i, st in enumerate(job.get("steps") or []):
            u = st.get("uses")
            if u is not None and not u.startswith("./") and not SHA.match(u):
                out.append(("W3", "job %s step %d: `uses: %s` is not pinned to a 40-hex commit SHA" % (name, i + 1, u)))
    if re.search(r"-latest\b", re.sub(r"#[^\n]*", "", text)):
        out.append(("W10", "a `-latest` runner label"))
    for m in re.finditer(r"^\s*-?\s*uses:\s*(\S+)[^\n]*$", text, re.M):
        if not m.group(1).startswith("./") and "#" not in m.group(0):
            out.append(("W3", "`%s` has no version comment" % m.group(0).strip()))
    if "pull_request_target" in text:
        out.append(("W4", "pull_request_target is used"))
    if re.search(r"secrets\.", text):
        out.append(("W4", "the workflow uses a secret"))
    if not re.search(r"sed -n 's/\^cancho \*= \*\"[^']*' cancho\.toml", text):
        out.append(("W5", "no step reads the compiler revision from cancho.toml with sed"))
    for m in re.finditer(r"(?<![0-9a-f])[0-9a-f]{40}(?![0-9a-f])", text):
        if not re.search(r"@" + m.group(0), text):
            out.append(("W5", "a 40-hex literal that is not an action pin: %s" % m.group(0)))
    for s in sorted(set(re.findall(r"(?<![\w./-])((?:scripts|ci)/[\w./-]+)", text))):
        s = s.rstrip(".")
        if not (root / s).exists():
            out.append(("W6", "the workflow names %s, which is not a file of the repository" % s))
    on = d.get("on") or {}
    for t in ("push", "pull_request", "workflow_dispatch"):
        if t not in on:
            out.append(("W7", "the trigger %s is missing" % t))
    fg = jobs.get("full-gate")
    if fg is None:
        out.append(("W7", "no full-gate job"))
    else:
        cond = str(fg.get("if", ""))
        if "workflow_dispatch" not in cond or "schedule" not in cond or "push" in cond or "pull_request" in cond:
            out.append(("W7", "full-gate must run only for workflow_dispatch and schedule (its `if` is %r)" % cond))
    bt = jobs.get("build-test") or {}
    inc = ((bt.get("strategy") or {}).get("matrix") or {}).get("include") or []
    oses = [str(i.get("os", "")) for i in inc]
    if not any(o.startswith("ubuntu") and not o.endswith("-arm") for o in oses) or not any(o.startswith("ubuntu") and o.endswith("-arm") for o in oses) or not any(o.startswith("macos") for o in oses):
        out.append(("W8", "build-test's matrix lacks a Linux x86-64, a Linux aarch64 or a macOS runner (%s)" % oses))
    return out


def mutants(text):
    def rep(a, b, n=1):
        assert a in text, a
        return text.replace(a, b, n)
    return [
        ("top-level permissions removed", rep("permissions:\n  contents: read\n\nconcurrency", "concurrency"), "W1"),
        ("top-level permissions write-all", rep("permissions:\n  contents: read\n\nconcurrency", "permissions: write-all\n\nconcurrency"), "W1"),
        ("top-level contents: write", rep("permissions:\n  contents: read\n\nconcurrency", "permissions:\n  contents: write\n\nconcurrency"), "W1"),
        ("a job writes", rep("    timeout-minutes: 30\n", "    timeout-minutes: 30\n    permissions:\n      contents: write\n"), "W2"),
        ("an action pinned by tag", rep("actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1 # v7.0.1", "actions/checkout@v4"), "W3"),
        ("an action pinned by a short SHA", rep("actions/cache@55cc8345863c7cc4c66a329aec7e433d2d1c52a9", "actions/cache@55cc834"), "W3"),
        ("an action pinned by a branch", rep("actions/setup-python@5fda3b95a4ea91299a34e894583c3862153e4b97", "actions/setup-python@main"), "W3"),
        ("an action pin with no version comment", rep("actions/upload-artifact@cf430e030ddbb5b0abf93d22962f4752f3646cd9 # v7.0.2", "actions/upload-artifact@cf430e030ddbb5b0abf93d22962f4752f3646cd9"), "W3"),
        ("pull_request_target", rep("  pull_request:\n", "  pull_request_target:\n"), "W4"),
        ("a secret", rep("persist-credentials: false", "token: ${{ secrets.DEPLOY_KEY }}"), "W4"),
        ("the compiler revision as a literal", rep("rev=$(sed -n 's/^cancho *= *\"\\([0-9a-f]*\\)\".*/\\1/p' cancho.toml)", "rev=de4e52fd81809e6400e2c8b2f485d0b75412449c", 1), "W5"),
        ("a revision literal in the environment", rep("concurrency:", "env:\n  CANCHO_REV: de4e52fd81809e6400e2c8b2f485d0b75412449c\n\nconcurrency:"), "W5"),
        ("a script that does not exist", rep("python3 scripts/check_workflow.py\n", "python3 scripts/check_nothing.py\n"), "W6"),
        ("a requirements file that does not exist", rep("ci/requirements.txt", "ci/requirements-none.txt", 1), "W6"),
        ("workflow_dispatch removed", rep("  workflow_dispatch:\n", ""), "W7"),
        ("pull_request removed", rep("  pull_request:\n", ""), "W7"),
        ("full-gate runs on every push", rep("if: github.event_name == 'workflow_dispatch' || github.event_name == 'schedule'", "if: github.event_name == 'workflow_dispatch' || github.event_name == 'schedule' || github.event_name == 'push'"), "W7"),
        ("no macOS in the matrix", rep("          - {os: macos-26, label: darwin-arm64}\n", ""), "W8"),
        ("no Linux aarch64 in the matrix", rep("          - {os: ubuntu-24.04-arm, label: linux-aarch64-ubuntu24}\n", ""), "W8"),
        ("a job without a timeout", rep("    timeout-minutes: 30\n", "", 1), "W9"),
        ("a -latest runner", rep("    runs-on: ubuntu-24.04\n    timeout-minutes: 30", "    runs-on: ubuntu-latest\n    timeout-minutes: 30"), "W10"),
        ("a runner without a version", rep("          - {os: macos-26, label: darwin-arm64}", "          - {os: macos-latest, label: darwin-arm64}"), "W10"),
        ("not YAML", text + "\n  - : [\n", "W0"),
    ]


def main(argv):
    mut = "--mutants" in argv
    args = [a for a in argv if not a.startswith("--")]
    path = pathlib.Path(args[0]) if args else DEFAULT
    if not path.exists():
        print("check_workflow: no %s" % path, file=sys.stderr)
        return 2
    text = path.read_text()
    found = check(text)
    if found:
        print("check_workflow: FAIL")
        for r, m in found:
            print("  %s %s" % (r, m))
        return 1
    print("check_workflow: ok (%s)" % path)
    if not mut:
        return 0
    bad = 0
    ms = mutants(text)
    for label, mt, rule in ms:
        got = {r for r, _ in check(mt)}
        ok = rule in got
        print("%-4s %s (%s)" % ("ok" if ok else "FAIL", label, rule))
        bad += 0 if ok else 1
    print("check_workflow --mutants: %d mutants, %d survived" % (len(ms), bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
