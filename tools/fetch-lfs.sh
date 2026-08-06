#!/usr/bin/env bash
#
# Fetches this checkout's Git LFS content, degrading gracefully instead of
# failing the job outright when GitHub's per-account LFS bandwidth/storage
# budget is exhausted for the current billing period.
#
# GitHub bills Git LFS bandwidth per ACCOUNT, not per job or per checkout.
# Every job in a CI run that fetches LFS content draws from the SAME shared
# budget, concurrently -- so once it runs dry mid-run, some matrix legs
# succeed (scheduled/fetched before the budget hit zero) and others fail
# with the identical server response. From the Actions UI that looks like a
# scattering of red jobs across every OS/leg with no code-level pattern to
# them, which is exactly what an exhausted-budget incident looks like: an
# account/billing condition, not something the PR under review broke.
# Failing CI red for it is actively misleading.
#
# The calling job is expected to have checked out with `lfs: false` (a fast,
# LFS-free clone -- see ci.yml) and to have already restored whatever it can
# from the `.git/lfs` object cache, so this script's own `git lfs pull` only
# needs to fetch what is genuinely new. This script performs that fetch and
# reports success/failure via `available` in $GITHUB_OUTPUT; the caller
# decides, via that output, whether to run the step(s) that need real
# fixture content -- see the `if: steps.lfs.outputs.available != 'false'`
# guards next to it in ci.yml.
#
# A genuine LFS outage or misconfiguration -- anything OTHER than the
# specific "exceeded its LFS budget" server response -- is NOT softened:
# this script only recognizes that one message. Every other failure exits
# non-zero and fails the job exactly as an unconditional `git lfs pull`
# always did.
#
# Usage: bash tools/fetch-lfs.sh
# Requires: $GITHUB_OUTPUT set (i.e. running inside a GitHub Actions job).

set -uo pipefail

git lfs install --local >/dev/null

echo "fetch-lfs: pulling LFS content..."
fetch_output=$(git lfs pull 2>&1)
status=$?

if [ "$status" -eq 0 ]; then
    echo "fetch-lfs: OK."
    echo "available=true" >>"$GITHUB_OUTPUT"
    exit 0
fi

if printf '%s\n' "$fetch_output" | grep -qi "exceeded its LFS budget"; then
    printf '%s\n' "$fetch_output"
    echo "::warning::This repository's Git LFS budget is exhausted for the current billing period, so this job's LFS-dependent step(s) are being SKIPPED rather than failed. This is an account/billing condition (GitHub organization/account Settings -> Billing -> Git LFS Data), not something this change broke -- re-run once the budget is restored to get real coverage from this job."
    echo "available=false" >>"$GITHUB_OUTPUT"
    exit 0
fi

echo "::error::git lfs pull failed for a reason OTHER than an exhausted LFS budget -- this is a real failure, not the degraded-budget case this script exists to soften:"
printf '%s\n' "$fetch_output"
exit "$status"
