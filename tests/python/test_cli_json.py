"""The two CLIs' ``--json`` outputs are one shape.

Every verb both CLIs share that prints JSON is run through the Python CLI
(in-process) and the native binary on the same files; each output must parse
as strict JSON (no bare ``NaN``), and the two must agree -- the same keys, the
same strings, integers, booleans and nulls, and floats within a tolerance -- as
must the exit codes. The native binary is found through
``MESHIOPLUSPLUS_NATIVE_CLI`` or under ``build/``; without one, the parity half
is skipped and only the Python half's validity is checked.
"""

import contextlib
import glob
import io
import json
import math
import os
import subprocess

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus._regions import Region

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Keys only the Python CLI reports (documented in doc/cli.md#json-output).
PYTHON_ONLY = {"sets"}


def _native_cli():
    env = os.environ.get("MESHIOPLUSPLUS_NATIVE_CLI")
    if env and os.path.isfile(env):
        return env
    found = [
        p
        for p in glob.glob(os.path.join(ROOT, "build", "*", "meshioplusplus"))
        if os.path.isfile(p) and os.access(p, os.X_OK)
    ]
    return max(found, key=os.path.getmtime) if found else None


NATIVE = _native_cli()


def _python(args):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = mio._cli.main(args)
    return code or 0, out.getvalue()


def _native(args):
    r = subprocess.run([NATIVE, *args], capture_output=True, text=True)
    return r.returncode, r.stdout


def _strict(text):
    def reject(token):
        raise ValueError(f"non-standard JSON constant {token}")

    return json.loads(text, parse_constant=reject)


def _same(a, b, path="$"):
    if isinstance(a, dict) and isinstance(b, dict):
        ka = set(a) - PYTHON_ONLY
        kb = set(b) - PYTHON_ONLY
        assert (
            ka == kb
        ), f"{path}: keys differ: python-only {ka - kb}, native-only {kb - ka}"
        for k in ka:
            _same(a[k], b[k], f"{path}.{k}")
    elif isinstance(a, list) and isinstance(b, list):
        assert len(a) == len(b), f"{path}: lengths {len(a)} != {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            _same(x, y, f"{path}[{i}]")
    elif isinstance(a, bool) or isinstance(b, bool) or a is None or b is None:
        assert a == b, f"{path}: {a!r} != {b!r}"
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12), f"{path}: {a} != {b}"
    else:
        assert a == b, f"{path}: {a!r} != {b!r}"


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    d = tmp_path_factory.mktemp("cli_json")
    n = 3
    xs, ys = np.meshgrid(np.arange(n + 1.0), np.arange(n + 1.0))
    pts = np.column_stack([xs.ravel(), ys.ravel(), np.zeros(xs.size)])
    quads = np.array(
        [
            [
                j * (n + 1) + i,
                j * (n + 1) + i + 1,
                (j + 1) * (n + 1) + i + 1,
                (j + 1) * (n + 1) + i,
            ]
            for j in range(n)
            for i in range(n)
        ]
    )
    grid = mio.Mesh(
        pts,
        [("quad", quads)],
        point_data={"u": pts[:, 0] * 2.0},
        cell_data={"c": [np.arange(len(quads), dtype=float)]},
    )
    moved = mio.Mesh(pts + [0.0, 0.0, 0.25], [("quad", quads)])
    tagged = mio.Mesh(pts, [("quad", quads)])
    tagged.regions = [
        Region("left", "point", [0, 4, 8, 12]),
        Region("right", "point", [3, 7, 11, 15]),
        Region("a", "cell", [0, 1, 2]),
        Region("b", "cell", [2, 3]),
    ]
    parts = mio.Mesh(
        np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1]],
            dtype=float,
        ),
        [
            (
                "tetra",
                np.array([[0, 1, 2, 3], [0, 2, 1, 4]], dtype=np.int64),
            )
        ],
        regions=[
            Region("upper", "cell", [0], dim=3),
            Region("lower", "cell", [1], dim=3),
            Region("slave", "point", [0]),
            Region("cut", "side", [[0, 3]]),
        ],
    )
    paths = {
        "grid": str(d / "grid.vtu"),
        "moved": str(d / "moved.vtu"),
        "tagged": str(d / "tagged.inp"),
        "parts": str(d / "parts.vtu"),
    }
    # a 3x3x3 hexahedral block, for agglomerate
    xs, ys, zs = np.meshgrid(*(np.arange(n + 1.0),) * 3, indexing="ij")
    hpts = np.column_stack([xs.ravel(), ys.ravel(), zs.ravel()])

    def vid(i, j, k):
        return (i * (n + 1) + j) * (n + 1) + k

    hexes = np.array(
        [
            [
                vid(i, j, k),
                vid(i + 1, j, k),
                vid(i + 1, j + 1, k),
                vid(i, j + 1, k),
                vid(i, j, k + 1),
                vid(i + 1, j, k + 1),
                vid(i + 1, j + 1, k + 1),
                vid(i, j + 1, k + 1),
            ]
            for i in range(n)
            for j in range(n)
            for k in range(n)
        ]
    )
    paths["hexes"] = str(d / "hexes.vtu")
    mio.write(paths["hexes"], mio.Mesh(hpts, [("hexahedron", hexes)]))
    # a two-step sequence (times 0 and 1 by index), for resample
    for k in range(2):
        mio.write(
            str(d / f"seq_{k}.vtu"),
            mio.Mesh(pts, [("quad", quads)], point_data={"u": pts[:, 0] + 10.0 * k}),
        )
    paths["seq"] = str(d / "seq_*.vtu")
    mio.write(paths["grid"], grid)
    mio.write(paths["moved"], moved)
    mio.write(paths["tagged"], tagged)
    mio.write(paths["parts"], parts)
    paths["dir"] = str(d)
    return paths


CASES = [
    ("info", lambda f: ["info", f["grid"], "--json"]),
    ("info_fast", lambda f: ["info", f["grid"], "--fast", "--json"]),
    ("info_regions", lambda f: ["info", f["tagged"], "--json"]),
    ("quality", lambda f: ["quality", f["grid"], "--json"]),
    (
        "check_pass",
        lambda f: ["check", f["grid"], "--require", "aspect_ratio <= 1.5", "--json"],
    ),
    (
        "check_fail",
        lambda f: ["check", f["grid"], "-r", "aspect_ratio <= 0.5 @ 10%", "--json"],
    ),
    ("stats", lambda f: ["stats", f["grid"], "--json"]),
    ("diff_same", lambda f: ["diff", f["grid"], f["grid"], "--json"]),
    ("diff_moved", lambda f: ["diff", f["grid"], f["moved"], "--json"]),
    ("regions", lambda f: ["regions", f["tagged"], "--json"]),
    ("data_info", lambda f: ["data", "info", f["grid"], "--json"]),
    (
        "hausdorff",
        lambda f: ["hausdorff", f["grid"], f["moved"], "--max", "0.1", "--json"],
    ),
    (
        "periodic",
        lambda f: [
            "periodic",
            f["tagged"],
            "--slave",
            "left",
            "--master",
            "right",
            "--translate",
            "3,0,0",
            "--json",
        ],
    ),
    (
        "contact_pairs",
        lambda f: [
            "contact-pairs",
            f["parts"],
            f["parts"],
            "--slave-region",
            "slave",
            "--master-region",
            "lower",
            "--tolerance",
            "10",
            "--json",
        ],
    ),
]

# Verbs that write a file: each CLI gets its own output path.
WRITING = [
    (
        "feature_edges",
        lambda f, out: ["feature-edges", f["grid"], out + ".vtu", "--json"],
    ),
    (
        "region_adjacency",
        lambda f, out: ["region-adjacency", f["tagged"], out + ".vtu", "--json"],
    ),
    (
        "find_interface",
        lambda f, out: [
            "find-interface",
            f["parts"],
            out + ".vtu",
            "--region-a",
            "upper",
            "--region-b",
            "lower",
            "--json",
        ],
    ),
    (
        "split_interface",
        lambda f, out: [
            "split-interface",
            f["parts"],
            out + ".vtu",
            "--side-region",
            "cut",
            "--add-cohesive",
            "--json",
        ],
    ),
    ("convert", lambda f, out: ["convert", f["grid"], out + ".vtk", "--json"]),
    (
        "agglomerate",
        lambda f, out: [
            "agglomerate",
            f["hexes"],
            out + ".vtu",
            "--merge-coplanar-faces",
            "--min-sphericity",
            "0.3",
            "--json",
        ],
    ),
    (
        "resample",
        lambda f, out: [
            "resample",
            f["seq"],
            out + "_{index}.vtu",
            "--time-from",
            "index",
            "--times",
            "0:1:0.25",
            "--json",
        ],
    ),
    (
        "resample_times_from",
        lambda f, out: [
            "resample",
            f["seq"],
            out + "_{index}.vtu",
            "--time-from",
            "index",
            "--times-from",
            f["seq"],
            "--json",
        ],
    ),
    (
        "resample_nearest",
        lambda f, out: [
            "resample",
            f["seq"],
            out + "_{index}.vtu",
            "--time-from",
            "index",
            "--times",
            "0.2,0.9",
            "--method",
            "nearest",
            "--json",
        ],
    ),
    (
        "regions_edit",
        lambda f, out: [
            "regions",
            f["tagged"],
            out + ".inp",
            "--union",
            "ab=a,b",
            "--json",
        ],
    ),
]


def _run_both(args_py, args_native):
    code_py, text_py = _python(args_py)
    doc_py = _strict(text_py)
    if NATIVE is None:
        return code_py, doc_py, None, None
    code_n, text_n = _native(args_native)
    return code_py, doc_py, code_n, _strict(text_n)


@pytest.mark.parametrize("name, build", CASES, ids=[c[0] for c in CASES])
def test_json_output_agrees(files, name, build):
    code_py, doc_py, code_n, doc_n = _run_both(build(files), build(files))
    if NATIVE is None:
        pytest.skip("no native CLI binary built; checked the Python output only")
    assert code_py == code_n
    _same(doc_py, doc_n)


@pytest.mark.parametrize("name, build", WRITING, ids=[w[0] for w in WRITING])
def test_json_output_of_writing_verbs_agrees(files, name, build):
    base = os.path.join(files["dir"], name)
    args_py = build(files, base + "_py")
    args_n = build(files, base + "_native")
    code_py, doc_py, code_n, doc_n = _run_both(args_py, args_n)
    if NATIVE is None:
        pytest.skip("no native CLI binary built; checked the Python output only")
    assert code_py == code_n
    for doc in (doc_py, doc_n):
        if isinstance(doc, dict) and "output" in doc:
            doc["output"] = "<out>"
    _same(doc_py, doc_n)


def test_json_is_strict_even_with_nan(files, capsys):
    # a metric that is NaN everywhere prints null, never a bare NaN
    code, text = _python(["quality", files["grid"], "--json"])
    doc = _strict(text)
    assert code == 0
    # min_angle is not defined for quads: count 0, and min/max/mean are null
    assert doc["metrics"]["quality:min_angle"]["count"] == 0
    assert doc["metrics"]["quality:min_angle"]["min"] is None


def test_check_exit_codes(files, tmp_path):
    assert _python(["check", files["grid"], "-r", "aspect_ratio <= 1.5"])[0] == 0
    assert _python(["check", files["grid"], "-r", "aspect_ratio <= 0.5"])[0] == 1
    assert _python(["check", files["grid"], "-r", "bogus >= 1"])[0] == 2
    assert _python(["check", str(tmp_path / "missing.vtu")])[0] == 2
    gate = tmp_path / "gate.txt"
    gate.write_text("# the part's gate\nskewness <= 0.1\naspect_ratio <= 1.5 @ 5%\n")
    assert _python(["check", files["grid"], "--gate", str(gate)])[0] == 0
    if NATIVE is not None:
        assert _native(["check", files["grid"], "-r", "aspect_ratio <= 0.5"])[0] == 1
        assert _native(["check", files["grid"], "-r", "bogus >= 1"])[0] == 2
        assert _native(["check", files["grid"], "--gate", str(gate)])[0] == 0
