#!/usr/bin/env python3

import argparse
import subprocess
import sys
from typing import List, Set


def run_git_command(args: List[str]) -> str:
    result = subprocess.run(
        ["git", *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        stderr = result.stderr.strip() or "unknown git error"
        raise RuntimeError(f"git {' '.join(args)} failed: {stderr}")

    return result.stdout


def get_local_branches() -> Set[str]:
    output = run_git_command(["for-each-ref", "--format=%(refname:short)", "refs/heads"])
    return {line.strip() for line in output.splitlines() if line.strip()}


def get_remote_branches(remote: str) -> Set[str]:
    output = run_git_command(["ls-remote", "--heads", remote])
    branch_names: Set[str] = set()

    for line in output.splitlines():
        parts = line.strip().split()
        if len(parts) != 2:
            continue

        ref = parts[1]
        if ref.startswith("refs/heads/"):
            branch_names.add(ref[len("refs/heads/") :])

    return branch_names


def get_current_branch() -> str:
    return run_git_command(["branch", "--show-current"]).strip()


def delete_branch(branch: str, force: bool) -> None:
    delete_flag = "-D" if force else "-d"
    run_git_command(["branch", delete_flag, branch])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Delete local branches that do not exist on a remote."
    )
    parser.add_argument(
        "--remote",
        default="origin",
        help="Remote name to compare against (default: origin)",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually delete stale branches (default is dry-run)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force delete branches (uses git branch -D)",
    )
    parser.add_argument(
        "--no-fetch",
        action="store_true",
        help="Do not run 'git fetch --prune <remote>' before checking",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        run_git_command(["rev-parse", "--is-inside-work-tree"])

        if not args.no_fetch:
            print(f"Fetching and pruning remote '{args.remote}'...")
            run_git_command(["fetch", "--prune", args.remote])

        local_branches = get_local_branches()
        remote_branches = get_remote_branches(args.remote)
        current_branch = get_current_branch()

        stale_branches = sorted(local_branches - remote_branches)

        if current_branch in stale_branches:
            stale_branches.remove(current_branch)
            print(
                f"Skipping current branch '{current_branch}' even though it is not on '{args.remote}'."
            )

        if not stale_branches:
            print("No stale local branches found.")
            return 0

        print("Stale local branches:")
        for branch in stale_branches:
            print(f"  - {branch}")

        if not args.apply:
            print("\nDry-run mode. Re-run with '--apply' to delete these branches.")
            return 0

        print("\nDeleting stale local branches...")
        for branch in stale_branches:
            delete_branch(branch, args.force)
            print(f"Deleted: {branch}")

        print("Done.")
        return 0

    except RuntimeError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
