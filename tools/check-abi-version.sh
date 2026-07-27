#!/usr/bin/env bash
#
# Keeps MESHIOPLUSPLUS_ABI_VERSION honest.
#
# A hand-bumped integer has exactly the failure mode of the prose rule it
# replaces: someone forgets. This is the gate that makes forgetting a red build.
#
#   Any installed header changed since the previous release tag
#   AND MESHIOPLUSPLUS_ABI_VERSION did not change
#   AND the release recorded no "reviewed as additive" entry
#   => FAIL
#
# It deliberately does NOT try to decide whether a header change was additive.
# That call needs a human -- editing the body of an existing inline function is a
# Tier B ODR break that looks textually identical to adding a brand new one -- so
# the tool demands the judgement be *recorded* in doc/abi_reviews.md rather than
# pretending to make it. `tests/cpp/test_abi_layout.cpp` is the complementary
# half: it catches Tier A layout changes mechanically, with no judgement at all.
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
    printf '%s\n' "$text" | sed -n 's/^#define[[:space:]]\+MESHIOPLUSPLUS_ABI_VERSION[[:space:]]\+\([0-9]\+\).*/\1/p' | head -1
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
# additive, and that review must be recorded against this ABI version.
if [ -f "$reviews" ] && grep -qE "^\|[[:space:]]*${abi_now}[[:space:]]*\|" "$reviews"; then
    if grep -qE "^\|[[:space:]]*${abi_now}[[:space:]]*\|.*\|" "$reviews"; then
        echo "check-abi-version: OK (headers changed, ABI held at ${abi_now}, and"
        echo "                      ${reviews} records the additive review)."
        exit 0
    fi
fi

cat >&2 <<EOF

check-abi-version: FAILED

  Installed headers changed between ${base} and HEAD, but
  MESHIOPLUSPLUS_ABI_VERSION is still ${abi_now} and ${reviews} records no
  review for it.

  Read doc/abi.md and decide which this was:

    Tier A (layout) or Tier B (an existing inline function/template body,
    or a default argument, changed)
        -> bump MESHIOPLUSPLUS_ABI_VERSION in ${abi_header}
           and update tests/cpp/test_abi_layout.cpp.

    Tier C (purely additive: a NEW inline function, constexpr, type,
    declaration or header, or an appended enumerator)
        -> add a row to ${reviews} for ABI ${abi_now} naming the headers
           and why the change cannot affect an already-compiled consumer.

  Not a rubber stamp: getting this wrong ships silent memory corruption to
  every C++ consumer that trusts the ABI number instead of the version pin.
EOF
exit 1
