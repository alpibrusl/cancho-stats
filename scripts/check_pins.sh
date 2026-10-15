#!/bin/sh
# G6 (docs/numerics.md section 2.1): the dependencies of cancho.toml are pinned stores, none copied. Ids:
#   PN1 a store of the pinned set has no block, or a block names a store that is not in the set
#   PN2 a label, or a (git, path), appears twice
#   PN3 a rev is not 40 hex digits, a repository's blocks carry different revs, or the rev is not the one the document names
#   PN4 a git is not the repository's GitHub URL
#   PN5 a path is not .cancho-vcs/<module>
#   PN6 a .cho file under lib, tools or tests declares a pinned module (a copy), or a sources entry of cancho.toml leaves the tree
#   PN7 the pinned set of the document (between the pins:begin and pins:end markers) differs from the set below
#   PN8 build/deps exists and lacks an installed module of the set
#   PN9 CHECK_REMOTE=1: the rev has the store on GitHub (skipped offline)
# TOML and DOC override the files read (the mutants). Exit 1 with the id and the reason on the first failure.
set -e
cd "$(dirname "$0")/.."
TOML=${TOML:-cancho.toml}
DOC=${DOC:-docs/numerics.md}
fail() { echo "check_pins: $1: $2" >&2; exit 1; }
[ -f "$TOML" ] || fail PN1 "no $TOML"
[ -f "$DOC" ] || fail PN7 "no $DOC"

# repository|url|doc-key
REPOS="cancho-tools|https://github.com/alpibrusl/cancho-tools
cancho-table|https://github.com/alpibrusl/cancho-table
cancho-ml|https://github.com/alpibrusl/cancho-ml"

SET=$(cat <<'EOF'
cancho-tools|cli|toolbox.cli
cancho-tools|describe|toolbox.describe
cancho-tools|fail|toolbox.fail
cancho-tools|limit|toolbox.limit
cancho-tools|lines|toolbox.lines
cancho-tools|out|toolbox.out
cancho-tools|path|toolbox.path
cancho-tools|place|toolbox.place
cancho-tools|text|toolbox.text
cancho-table|dec|dec
cancho-table|facc|facc
cancho-table|flt|flt
cancho-table|pow5|pow5
cancho-table|reader|reader
cancho-ml|ml-refuse|ml.refuse
cancho-ml|ml-limits|ml.limits
cancho-ml|ml-rng|ml.rng
cancho-ml|ml-exact|ml.exact
cancho-ml|ml-fbits|ml.fbits
cancho-ml|ml-fsq|ml.fsq
cancho-ml|ml-matrix|ml.matrix
cancho-ml|linalg-order|linalg.order
cancho-ml|linalg-tensor|linalg.tensor
cancho-ml|linalg-dot|linalg.dot
cancho-ml|linalg-axpy|linalg.axpy
cancho-ml|linalg-gemv|linalg.gemv
cancho-ml|linalg-syrk|linalg.syrk
cancho-ml|linalg-gemm|linalg.gemm
cancho-ml|linalg-chol|linalg.chol
cancho-ml|linalg-par|linalg.par
cancho-ml|metrics-rules|metrics.rules
cancho-ml|metrics-ratio|metrics.ratio
cancho-ml|metrics-classify|metrics.classify
cancho-ml|metrics-sums|metrics.sums
cancho-ml|metrics-regress|metrics.regress
cancho-ml|metrics-rank|metrics.rank
cancho-ml|metrics-splits|metrics.splits
EOF
)

# PN7: the document's list
doclist=$(awk '/<!-- pins:begin -->/{f=1;next} /<!-- pins:end -->/{f=0} f' "$DOC" | sed -n 's/^| *\([^ |]*\) *| *\([^ |]*\) *| *\([^ |]*\) *|.*/\1|\2|\3/p')
[ -n "$doclist" ] || fail PN7 "$DOC has no pinned set between the pins:begin and pins:end markers"
a=$(echo "$SET" | sort); b=$(echo "$doclist" | sort)
[ "$a" = "$b" ] || fail PN7 "the pinned set of $DOC differs from the set of this script: $(echo "$a
$b" | sort | uniq -u | head -3 | tr '\n' ' ')"

docrev() { sed -n "s/^$1 rev: \`\([0-9a-f]*\)\`.*/\1/p" "$DOC" | head -1; }
# the compiler
ccwant=$(docrev cancho)
[ ${#ccwant} -eq 40 ] || fail PN3 "$DOC has no line 'cancho rev: \`<40 hex digits>\`'"
cc=$(sed -n 's/^cancho *= *"\([0-9a-f]*\)".*/\1/p' "$TOML" | head -1)
[ "$cc" = "$ccwant" ] || fail PN3 "[package] cancho is '$cc', not the documented $ccwant"

# blocks of the toml: label|git|rev|path
blocks=$(awk '
  /^\[dependencies\./ { if (lab!="") print lab "|" git "|" rev "|" path; lab=$0; sub(/^\[dependencies\./,"",lab); sub(/\]$/,"",lab); git=""; rev=""; path=""; next }
  /^\[/ { if (lab!="") print lab "|" git "|" rev "|" path; lab=""; next }
  lab!="" && /^git *=/ { s=$0; sub(/^git *= *"/,"",s); sub(/".*$/,"",s); git=s }
  lab!="" && /^rev *=/ { s=$0; sub(/^rev *= *"/,"",s); sub(/".*$/,"",s); rev=s }
  lab!="" && /^path *=/ { s=$0; sub(/^path *= *"/,"",s); sub(/".*$/,"",s); path=s }
  END { if (lab!="") print lab "|" git "|" rev "|" path }
' "$TOML")

# PN2: duplicates
dup=$(echo "$blocks" | cut -d'|' -f1 | sort | uniq -d | head -1)
[ -z "$dup" ] || fail PN2 "[dependencies.$dup] is listed twice"
dup=$(echo "$blocks" | awk -F'|' '{print $2 "|" $4}' | sort | uniq -d | head -1)
[ -z "$dup" ] || fail PN2 "the same git and path appear twice: $dup"

# PN1: every block is in the set, every store of the set has a block
echo "$blocks" | while IFS='|' read -r lab git rev path; do
  [ -n "$lab" ] || continue
  echo "$SET" | cut -d'|' -f2 | grep -qx "$lab" || fail PN1 "[dependencies.$lab] is not in the pinned set (an unlisted dependency)"
done
echo "$SET" | while IFS='|' read -r repo lab mod; do
  echo "$blocks" | cut -d'|' -f1 | grep -qx "$lab" || fail PN1 "the pinned set has $lab ($mod), $TOML has no [dependencies.$lab]"
done

echo "$blocks" | while IFS='|' read -r lab git rev path; do
  [ -n "$lab" ] || continue
  entry=$(echo "$SET" | awk -F'|' -v l="$lab" '$2==l{print $1 "|" $3}')
  repo=${entry%%|*}; mod=${entry#*|}
  url=$(echo "$REPOS" | awk -F'|' -v r="$repo" '$1==r{print $2}')
  # PN3
  case "$rev" in *[!0-9a-f]*|"") fail PN3 "[dependencies.$lab]: rev '$rev' is not a full commit hash" ;; esac
  [ ${#rev} -eq 40 ] || fail PN3 "[dependencies.$lab]: rev is not a full commit hash (${#rev} digits)"
  want=$(docrev "$repo")
  [ ${#want} -eq 40 ] || fail PN3 "$DOC has no line '$repo rev: \`<40 hex digits>\`'"
  [ "$rev" = "$want" ] || fail PN3 "[dependencies.$lab]: rev $rev is not the documented $want ($repo)"
  # PN4
  [ "$git" = "$url" ] || fail PN4 "[dependencies.$lab]: git is '$git', not $url"
  # PN5
  [ "$path" = ".cancho-vcs/$mod" ] || fail PN5 "[dependencies.$lab]: path is '$path', not .cancho-vcs/$mod"
done

# one rev per repository (PN3)
for repo in cancho-tools cancho-table cancho-ml; do
  revs=$(echo "$blocks" | while IFS='|' read -r lab git rev path; do
    r=$(echo "$SET" | awk -F'|' -v l="$lab" '$2==l{print $1}'); [ "$r" = "$repo" ] && echo "$rev"; done | sort -u | wc -l | tr -d ' ')
  [ "$revs" -le 1 ] || fail PN3 "the blocks of $repo carry $revs different revs"
done

# PN6
if grep -q '^sources *=.*\.\./' "$TOML"; then fail PN6 "a sources entry of $TOML reaches outside this tree"; fi
mods=$(echo "$SET" | cut -d'|' -f3)
for f in $(find lib tools tests -name '*.cho' 2>/dev/null); do
  for m in $mods; do
    if grep -q "^module $m;" "$f"; then fail PN6 "$f declares module $m, which is a pinned store's"; fi
  done
done

# PN8
if [ -d build/deps ] && ls build/deps/*.cho >/dev/null 2>&1; then
  for m in $mods; do
    grep -qs "^module $m;" build/deps/*.cho || fail PN8 "build/deps has no installed module $m (run cancho install)"
  done
fi

# PN9
if [ -n "$CHECK_REMOTE" ]; then
  if command -v gh >/dev/null 2>&1; then
    echo "$SET" | while IFS='|' read -r repo lab mod; do
      rev=$(echo "$blocks" | awk -F'|' -v l="$lab" '$1==l{print $3}')
      gh api "repos/alpibrusl/$repo/contents/.cancho-vcs/$mod?ref=$rev" >/dev/null 2>&1 || fail PN9 "$repo $rev has no .cancho-vcs/$mod on GitHub"
    done
    echo "check_pins: remote ok"
  else
    echo "check_pins: remote check skipped (no gh)"
  fi
fi
n=$(echo "$SET" | wc -l | tr -d ' ')
echo "check_pins: ok ($n dependencies, none copied)"
