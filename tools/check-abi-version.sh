#!/usr/bin/env bash
#
# Keeps MESHIOPLUSPLUS_ABI_VERSION honest.
#
# A hand-bumped integer has exactly the failure mode of the prose rule it
# replaces: someone forgets. This is the gate that makes forgetting a red build.
#
#   Any installed header changed since the previous release tag
#   AND MESHIOPLUSPLUS_ABI_VERSION did not change
#   AND doc/abi_reviews.md has no row for THIS release naming THOSE headers
#   => FAIL
#
# It deliberately does NOT try to decide whether a header change was additive.
# That call needs a human -- editing the body of an existing inline function is a
# Tier B ODR break that looks textually identical to adding a brand new one -- so
# the tool demands the judgement be *recorded* in doc/abi_reviews.md rather than
# pretending to make it. `tests/cpp/test_abi_layout.cpp` is the complementary
# half: it catches Tier A layout changes mechanically, with no judgement at all.
#
# ### What the review lookup keys on, and why it is not just the ABI number
#
# Through v9.4.0 this matched any row whose first column was the current ABI
# number. That passes VACUOUSLY the moment one such row exists: doc/abi_reviews.md
# has carried two ABI-3 rows since v9.4.0, so every later header change matched
# one of them and exited 0 announcing "records the additive review" -- naming a
# review of an entirely different release. A v9.5.0 that edited an existing
# inline function body, a genuine Tier B break, would have sailed through: the
# layout snapshot cannot see an inline-body edit, and this was the only other
# gate.
#
# Keying on the release version alone does not fix it either. A header change
# that does not ALSO bump the version leaves $version_now at the previous
# release, which by construction already has a row. So the row must additionally
# account for the headers actually touched -- which is what the "headers changed"
# column of doc/abi_reviews.md is for, and why that column is load-bearing rather
# than decorative.
#
# See doc/abi.md for the tiers.
#
# Usage:
#   tools/check-abi-version.sh [BASE_REF]
#
# BASE_REF defaults to the most recent v* tag that is an ancestor of HEAD.
# CI must check out with tags available (actions/checkout fetch-depth: 0).

set -euo pipefail

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"

abi_header="src/cpp/include/meshioplusplus/abi_version.hpp"
reviews="doc/abi_reviews.md"
header_dir="src/cpp/include/meshioplusplus"

read_abi() {  # read_abi <ref|WORKTREE>
    local ref="$1" text
    if [ "$ref" = "WORKTREE" ]; then
        text=$(cat "$abi_header")
    else
        text=$(git show "$ref:$abi_header" 2>/dev/null || true)
    fi
    # `sed -nE`, not `sed -n` with `\+`: `\+` is a GNU extension that BSD sed
    # (macOS) reads as a literal plus, so the pattern silently matched nothing
    # there and the script died with "cannot parse MESHIOPLUSPLUS_ABI_VERSION".
    # -E/ERE is understood by GNU, BSD and busybox alike.
    printf '%s\n' "$text" | sed -nE 's/^#define[[:space:]]+MESHIOPLUSPLUS_ABI_VERSION[[:space:]]+([0-9]+).*/\1/p' | head -1
}

base=${1:-}
if [ -z "$base" ]; then
    base=$(git describe --tags --abbrev=0 --match 'v*' HEAD 2>/dev/null || true)
    # On the release commit ITSELF `git describe` returns that very tag, which
    # would compare the release against itself and pass vacuously. Step back one
    # -- but only when the working tree is clean. With uncommitted work on top of
    # a tagged HEAD (a release being prepared), that tag IS the right baseline,
    # and stepping back would silently measure against the wrong release.
    if [ -n "$base" ] && [ "$(git rev-parse "$base^{commit}")" = "$(git rev-parse HEAD)" ] &&
       [ -z "$(git status --porcelain -- "$header_dir" "$abi_header")" ]; then
        base=$(git describe --tags --abbrev=0 --match 'v*' HEAD~1 2>/dev/null || true)
    fi
fi

if [ -z "$base" ]; then
    echo "check-abi-version: no v* tag reachable from HEAD; nothing to compare against."
    echo "  (CI needs actions/checkout with fetch-depth: 0 for this gate to mean anything.)"
    exit 0
fi

echo "check-abi-version: comparing installed headers ${base}..HEAD"

# `$base --` and not `$base..HEAD`: this compares the tag against the WORKING
# TREE, so uncommitted work is measured too. In CI the tree is clean and the two
# forms agree; locally, `..HEAD` would quietly ignore the very change being made,
# which is the one moment the gate is supposed to speak up.
changed=$(git diff --name-only "$base" -- "$header_dir" || true)
if [ -z "$changed" ]; then
    echo "check-abi-version: no installed header changed. Nothing to check."
    exit 0
fi

abi_before=$(read_abi "$base")
abi_now=$(read_abi WORKTREE)

if [ -z "$abi_now" ]; then
    echo "check-abi-version: ERROR: cannot parse MESHIOPLUSPLUS_ABI_VERSION from $abi_header" >&2
    exit 1
fi

echo "check-abi-version: ABI ${abi_before:-<absent>} -> ${abi_now}"
echo "check-abi-version: installed headers changed:"
printf '  %s\n' $changed

# abi_version.hpp did not exist before it was introduced; that release is by
# definition the baseline and has nothing to compare.
if [ -z "$abi_before" ]; then
    echo "check-abi-version: OK ($abi_header did not exist at $base; baseline release)."
    exit 0
fi

if [ "$abi_before" != "$abi_now" ]; then
    echo "check-abi-version: OK (ABI version was bumped)."
    exit 0
fi

# Unchanged ABI + changed headers => the change must have been reviewed as
# additive, and that review must be recorded against THIS release and name THESE
# headers. See the header comment for why neither key alone is enough.
version_now=$(sed -nE 's/^[[:space:]]*VERSION[[:space:]]+([0-9][0-9.]*)[[:space:]]*$/\1/p' \
              CMakeLists.txt | head -1)
if [ -z "$version_now" ]; then
    echo "check-abi-version: ERROR: cannot parse the project VERSION from CMakeLists.txt" >&2
    exit 1
fi

# Rows for (this ABI, this release). Field-splitting on `|` with awk rather than
# one in-column regex: the release cell is free text, and matching a version
# inside it with grep needs an escaping dance that is easy to get subtly wrong.
# A release matches as a whole token, with or without the `v` prefix, so both
# "v9.4.1" and "9.4.1" read correctly and "9.4.10" does not match "9.4.1".
rows=$(awk -F'|' -v abi="$abi_now" -v ver="$version_now" '
    /^[[:space:]]*\|/ {
        field = $2
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", field)
        if (field != abi) next
        n = split($3, toks, /[^0-9A-Za-z.]+/)
        for (i = 1; i <= n; i++) {
            tok = toks[i]
            sub(/^v/, "", tok)
            if (tok == ver) { print; next }
        }
    }' "$reviews" 2>/dev/null || true)

if [ -n "$rows" ]; then
    # Every changed header must be accounted for by name. Paths are recorded
    # relative to the include root, the way doc/abi_reviews.md already writes
    # them (`formats/exodus.hpp`, not the full src/cpp/include/... path).
    missing=""
    for header in $changed; do
        rel=${header#"$header_dir/"}
        if ! printf '%s\n' "$rows" | grep -qF "$rel"; then
            missing="$missing $rel"
        fi
    done

    if [ -z "$missing" ]; then
        echo "check-abi-version: OK (headers changed, ABI held at ${abi_now}, and"
        echo "                      ${reviews} reviews them for v${version_now})."
        exit 0
    fi

    cat >&2 <<EOF

check-abi-version: FAILED

  ${reviews} has a row for ABI ${abi_now} / v${version_now}, but it does not
  account for every installed header that changed since ${base}:
$(printf '    %s\n' $missing)

  Either name them in that row -- having checked each one is Tier C by
  doc/abi.md's criterion, i.e. that it cannot affect an already-compiled
  consumer -- or bump MESHIOPLUSPLUS_ABI_VERSION in ${abi_header} if any of them
  was in fact a Tier A layout change or a Tier B inline-body/default-argument
  change.

  Not a rubber stamp: getting this wrong ships silent memory corruption to
  every C++ consumer that trusts the ABI number instead of the version pin.
EOF
    exit 1
fi

cat >&2 <<EOF

check-abi-version: FAILED

  Installed headers changed between ${base} and HEAD, but
  MESHIOPLUSPLUS_ABI_VERSION is still ${abi_now} and ${reviews} records no
  review for ABI ${abi_now} / v${version_now}.

  (A row for ABI ${abi_now} against some *other* release does not count. The
  review has to be of this change.)

  Read doc/abi.md and decide which this was:

    Tier A (layout) or Tier B (an existing inline function/template body,
    or a default argument, changed)
        -> bump MESHIOPLUSPLUS_ABI_VERSION in ${abi_header}
           and update tests/cpp/test_abi_layout.cpp.

    Tier C (purely additive: a NEW inline function, constexpr, type,
    declaration or header, or an appended enumerator)
        -> add a row to ${reviews} for ABI ${abi_now} / v${version_now},
           naming every header above and why the change cannot affect an
           already-compiled consumer.

  Not a rubber stamp: getting this wrong ships silent memory corruption to
  every C++ consumer that trusts the ABI number instead of the version pin.
EOF
exit 1
