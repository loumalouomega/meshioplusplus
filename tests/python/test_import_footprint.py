"""What ``import meshioplusplus`` pulls in (roadmap §3, "Boundaries and startup").

Deterministic, unlike a timing: each module below was measured on the import
path once and removed from it, and this keeps it out. Run in a subprocess so
the test sees a fresh interpreter, not whatever the suite imported already.
"""

import importlib.util
import json
import subprocess
import sys

# Modules the package must not import on its own:
#   rich -- ~17 ms; only the CLI and the three stderr helpers in _common use it.
#   xml.sax.saxutils and what it drags in (urllib.request, http.client, ssl)
#     -- ~12 ms, for one attribute-quoting function in _pvtk_index.
#   meshioplusplus._cli -- ~10 ms; only the command line needs it, and the
#     package loads it on first attribute access (see test below).
_ABSENT = [
    "rich",
    "xml.sax.saxutils",
    "urllib.request",
    "http.client",
    "ssl",
    "meshioplusplus._cli",
]


def _imported_after_import(names):
    """``{name: imported}`` after ``import meshioplusplus`` in a fresh
    interpreter; a module some ``.pth`` hook (an editable install's finder)
    had already imported before the package counts as ``None``, unknown."""
    code = (
        "import json, sys\n"
        f"before = {{n: n in sys.modules for n in {names!r}}}\n"
        "import meshioplusplus\n"
        f"print(json.dumps({{n: None if before[n] else n in sys.modules for n in {names!r}}}))\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], check=True, capture_output=True, text=True
    )
    return json.loads(out.stdout.strip().splitlines()[-1])


def test_heavy_modules_are_not_imported():
    present = _imported_after_import(_ABSENT)
    assert not any(present.values()), {n: p for n, p in present.items() if p}


def test_version_does_not_need_importlib_metadata():
    # An installed build carries a generated _version.py; a bare source tree
    # on sys.path falls back to importlib.metadata, which is fine there.
    if importlib.util.find_spec("meshioplusplus._version") is None:
        return
    present = _imported_after_import(["importlib.metadata"])
    assert present["importlib.metadata"] is not True


def test_cli_loads_on_first_attribute_access_and_through_star_import():
    code = (
        "import sys\n"
        "import meshioplusplus\n"
        "assert 'meshioplusplus._cli' not in sys.modules\n"
        "assert callable(meshioplusplus._cli.main)\n"
        "assert 'meshioplusplus._cli' in sys.modules\n"
        "ns = {}\n"
        "exec('from meshioplusplus import *', ns)\n"
        "assert ns['_cli'] is meshioplusplus._cli\n"
        "try:\n"
        "    meshioplusplus.no_such_attribute\n"
        "except AttributeError:\n"
        "    pass\n"
        "else:\n"
        "    raise SystemExit('missing attribute did not raise')\n"
    )
    subprocess.run([sys.executable, "-c", code], check=True, capture_output=True)


def test_extension_order_does_not_depend_on_the_cli_being_imported_first():
    # `_cli` used to be imported first and pulled its formats in ahead of the
    # alphabetical list, which is what put ansys, gmsh, freefem in that order.
    code = (
        "import json, meshioplusplus\n"
        "from meshioplusplus._helpers import extension_to_filetypes as m\n"
        "print(json.dumps(m['.msh']))\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], check=True, capture_output=True, text=True
    )
    assert json.loads(out.stdout.strip().splitlines()[-1]) == [
        "ansys",
        "gmsh",
        "freefem",
    ]
