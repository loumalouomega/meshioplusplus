"""Fail when the public Python API changed incompatibly without a ``**Breaking:**`` line.

Compares ``tools/api/python_api.json`` (see ``tools/gen_api_dump.py``) at a base
revision with the working tree. These count as breaking:

* a name leaving ``__all__``, or changing kind (function, class, constant);
* a parameter removed or renamed, a parameter that became required, a changed
  default, or a changed parameter kind;
* the same rules for the public methods of a class.

Pure additions, including a new parameter that has a default, are fine. When the
comparison finds something, ``CHANGELOG.md`` must have gained a line that contains
``**Breaking:**`` since the base::

    python tools/check_api_changes.py [BASE_REF]    # default: origin/master's merge base
"""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
DUMP = "tools/api/python_api.json"


def _params(signature):
    if not signature:
        return None
    return {p["name"]: p for p in signature["parameters"]}


def _compare_signature(where, old, new, out):
    po, pn = _params(old), _params(new)
    if po is None or pn is None:
        return  # unreadable on one side: nothing to compare
    for name, p in po.items():
        q = pn.get(name)
        if q is None:
            out.append(f"{where}: parameter '{name}' removed or renamed")
            continue
        if p["kind"] != q["kind"]:
            out.append(
                f"{where}: parameter '{name}' changed kind {p['kind']} -> {q['kind']}"
            )
        if p.get("default") != q.get("default"):
            out.append(
                f"{where}: default of '{name}' changed "
                f"{p.get('default', '<required>')} -> {q.get('default', '<required>')}"
            )
    for name, q in pn.items():
        if (
            name not in po
            and "default" not in q
            and q["kind"] not in ("VAR_POSITIONAL", "VAR_KEYWORD")
        ):
            out.append(f"{where}: new required parameter '{name}'")


def compare(old, new):
    """The list of incompatible changes from dump ``old`` to dump ``new``."""
    out = []
    names_old, names_new = old.get("names", {}), new.get("names", {})
    for name, eo in sorted(names_old.items()):
        en = names_new.get(name)
        if en is None:
            out.append(f"{name}: removed from __all__")
            continue
        if "unavailable" in (eo["kind"], en["kind"]):
            continue
        if eo["kind"] != en["kind"]:
            out.append(f"{name}: changed from {eo['kind']} to {en['kind']}")
            continue
        _compare_signature(name, eo.get("signature"), en.get("signature"), out)
        for member, mo in sorted(eo.get("members", {}).items()):
            mn = en.get("members", {}).get(member)
            if mn is None:
                out.append(f"{name}.{member}: removed")
            elif mo.get("kind") == "method" and mn.get("kind") == "method":
                _compare_signature(
                    f"{name}.{member}", mo.get("signature"), mn.get("signature"), out
                )
    return out


def _git(*args):
    return subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True)


def _base(argv):
    if argv:
        return argv[0]
    for ref in ("origin/master", "master"):
        r = _git("merge-base", "HEAD", ref)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip()
    return None


def declared_breaking(base):
    """True when CHANGELOG.md gained a ``**Breaking:**`` line since ``base``."""
    r = _git("diff", "-U0", base, "--", "CHANGELOG.md")
    return any(
        line.startswith("+") and not line.startswith("+++") and "**Breaking:**" in line
        for line in r.stdout.splitlines()
    )


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    base = _base(argv)
    if base is None:
        print("check_api_changes: no base revision to compare with; nothing to check")
        return 0
    old_text = _git("show", f"{base}:{DUMP}")
    if old_text.returncode != 0:
        print(
            f"check_api_changes: {DUMP} does not exist at {base[:10]}; nothing to compare"
        )
        return 0
    new_path = REPO / DUMP
    if not new_path.is_file():
        print(
            f"check_api_changes: {DUMP} is missing from the working tree",
            file=sys.stderr,
        )
        return 1
    changes = compare(json.loads(old_text.stdout), json.loads(new_path.read_text()))
    if not changes:
        print("check_api_changes: no incompatible API change")
        return 0
    print("check_api_changes: incompatible changes since", base[:10])
    for c in changes:
        print("  -", c)
    if declared_breaking(base):
        print("check_api_changes: declared with a **Breaking:** line in CHANGELOG.md")
        return 0
    print(
        "check_api_changes: add a **Breaking:** entry to CHANGELOG.md for these",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
