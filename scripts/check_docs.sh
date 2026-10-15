#!/bin/sh
# G10 (docs/numerics.md section 7.4): no tracked file over 2,000 lines; the documents name only scripts, tools and libraries that exist; and they carry no changelog narration
# (the phrases below belong to a history, and the documents state what is true now). A phrase can appear only in a line that quotes it between backticks.
set -e
cd "$(dirname "$0")/.."
fail() { echo "check_docs: $1" >&2; exit 1; }
FILES=$(git ls-files 2>/dev/null || find . -type f -not -path './.git/*' -not -path './build/*' -not -path './.venv/*')
for f in $FILES; do
  [ -f "$f" ] || continue
  n=$(wc -l < "$f" | tr -d ' ')
  [ "$n" -le 2000 ] || fail "$f has $n lines (the limit is 2,000)"
done
# docs/design.md is the maintainer's design of the whole tool and names the files of other repositories; the checks below are for the documents of this milestone
for d in docs/numerics.md; do
  [ -f "$d" ] || continue
  for ref in $(grep -o '\(scripts\|tools\|lib\|ci\)/[A-Za-z0-9_./-]*[A-Za-z0-9_]' "$d" | sort -u); do
    case "$ref" in *'*'*) continue ;; esac
    [ -e "$ref" ] || fail "$d names $ref, which does not exist"
  done
  hit=$(grep -n -i -E 'was added|were added|has been added|have been added|used to be|previously|no longer|changelog|now supports|newly added|we changed|we fixed' "$d" | grep -v '`[^`]*\(was added\|previously\|no longer\|changelog\)[^`]*`' || true)
  [ -z "$hit" ] || fail "$d has changelog narration: $(echo "$hit" | head -2 | tr '\n' ' ')"
done
echo "check_docs: ok"
