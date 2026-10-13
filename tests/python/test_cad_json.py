"""Tests for Kratos CAD JSON (``*.cad.json``), the NURBS B-rep of Kratos's IGA.

The fixtures in ``meshes/cad_json/`` come from ``tools/gen_cad_json_fixtures.py``
(closed-form geometry, written with meshio++'s own writer). Both engines must
read the same mesh and side channel and write the same bytes; the geometry must
be exact where the NURBS is (circles, a cylinder). An optional check (skipped
without ``KratosMultiphysics`` and its IgaApplication) loads what meshio++
writes through Kratos's own ``CadJsonInput``. The C++ behaviour is pinned in
``tests/cpp/test_cad_json.cpp``.
"""

import json
import math
import pathlib

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _json_text
from meshioplusplus import _nurbs_eval as ne
from meshioplusplus._exceptions import ReadError, WriteError
from meshioplusplus.cad_json import _cad_json as ref

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")

MESHES = pathlib.Path(__file__).parent / "meshes" / "cad_json"
FIXTURES = sorted(p.name for p in MESHES.glob("*.json"))
CANONICAL = [n for n in FIXTURES if n.endswith(".cad.json")]


def _same_mesh(a, b):
    assert np.array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        assert np.array_equal(x.data, y.data)
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for k in a.cell_data:
        for x, y in zip(a.cell_data[k], b.cell_data[k]):
            assert np.array_equal(x, y), k


def _same_brep(a, b):
    """Two side channels are equal value for value (arrays by content)."""
    if isinstance(a, dict):
        assert sorted(a) == sorted(b)
        for k in a:
            _same_brep(a[k], b[k])
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            _same_brep(x, y)
    elif isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        assert np.array_equal(np.asarray(a), np.asarray(b))
    else:
        assert a == b


# ---------------------------------------------------------------------------
# The kernels
# ---------------------------------------------------------------------------


def test_json_text_writes_the_nlohmann_dump4_layout():
    value = _json_text.parse(
        '{"a": [1, 2.5, {}], "b": [], "c": {"d": "x\\ny"}, "e": true}', "t"
    )
    assert _json_text.dumps(value) == (
        '{\n    "a": [\n        1,\n        2.5,\n        {}\n    ],\n    "b": [],\n'
        '    "c": {\n        "d": "x\\ny"\n    },\n    "e": true\n}'
    )
    assert (
        _json_text.dumps(value, indent=0)
        == '{"a":[1,2.5,{}],"b":[],"c":{"d":"x\\ny"},"e":true}'
    )


@pytest.mark.parametrize(
    "value, text",
    [
        (0.1, "0.1"),
        (15.0, "15"),
        (-0.0, "-0"),
        (1e-5, "1e-05"),
        (1 / 3, "0.3333333333333333"),
        (2.0**60, "1.152921504606847e+18"),
        (math.inf, "null"),
    ],
)
def test_json_numbers_read_back_exactly(value, text):
    assert _json_text.format_number(value) == text
    if math.isfinite(value):
        assert float(text) == value


def test_json_comments_order_and_duplicates():
    v = _json_text.parse('// c\n{"b": 1, /* x */ "a": 2, "b": 3}', "t")
    assert list(v) == [("b", 1), ("a", 2), ("b", 3)]
    assert v.find("b") == 1


@pytest.mark.parametrize(
    "text", ["{", '{"a": NaN}', "[1,]", '{"a" 1}', "[" * 300 + "]" * 300, "1 2"]
)
def test_json_rejects_malformed_text(text):
    with pytest.raises(ReadError):
        _json_text.parse(text, "t")


def test_nurbs_kernel_partition_of_unity_and_exact_circle():
    knots = [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]
    for u in np.linspace(0, 1, 17):
        span = ne.find_span(2, knots, 4, float(u))
        assert abs(sum(ne.basis(span, float(u), 2, knots)) - 1.0) < 1e-15
    w = math.sqrt(0.5)
    cps = [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    for u in np.linspace(0, 1, 11):
        x, y, _ = ne.curve_point_in_span(
            2, [0, 0, 0, 1, 1, 1], cps, [1.0, w, 1.0], 2, float(u)
        )
        assert abs(math.hypot(x, y) - 1.0) < 1e-15


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


@needs_core
@pytest.mark.parametrize("name", FIXTURES)
@pytest.mark.parametrize("subdivisions", [0, 1, 3])
def test_engines_agree(name, subdivisions):
    a = mp.cad_json.read(MESHES / name, subdivisions=subdivisions)
    b = ref.read(MESHES / name, subdivisions=subdivisions)
    _same_mesh(a, b)
    _same_brep(a.info["brep"], b.info["brep"])


@pytest.mark.parametrize("name", FIXTURES)
def test_dispatch_and_sniffing(name):
    path = MESHES / name
    mesh = (
        mp.read(path)
        if name.endswith(".cad.json")
        else mp.read(path, file_format="cad_json")
    )
    assert "brep" in mesh.info
    from meshioplusplus._sniff import _sniff_format_py

    assert _sniff_format_py(path) == "cad_json"
    if _core is not None:
        assert _core.sniff_format(str(path)) == "cad_json"


def test_cylinder_is_exact():
    mesh = mp.cad_json.read(MESHES / "cylindrical_roof.cad.json", subdivisions=8)
    p = mesh.points
    assert np.abs(np.hypot(p[:, 0], p[:, 2]) - 25.0).max() < 1e-12
    face = mesh.info["brep"]["faces"][0]
    assert face["rational"] and face["reduced_knots"] and not face["trimmed"]
    # Kratos's reduced knots are held in the full form.
    assert len(mesh.info["brep"]["patches"][0]["knots"][1]) == 4


def test_trimmed_plate_drops_the_hole():
    mesh = mp.cad_json.read(MESHES / "plate_with_hole.cad.json", subdivisions=8)
    faces = [i for i, c in enumerate(mesh.cells) if c.type.startswith("quad")]
    assert faces
    for i in faces:
        centres = mesh.points[mesh.cells[i].data[:, :4]].mean(axis=1)
        assert (np.hypot(centres[:, 0] - 2, centres[:, 1] - 2) > 0.85).all()
        trimmed = mesh.cell_data["cad:trimmed"][i]
        # Only cells along the hole are crossed; the outer boundary is not.
        crossed = centres[trimmed == 1]
        assert (
            len(crossed)
            and (np.hypot(crossed[:, 0] - 2, crossed[:, 1] - 2) < 1.5).all()
        )
    # The hole's edges are exact circles.
    edge_pts = np.unique(
        np.concatenate(
            [
                c.data.ravel()
                for c, k in zip(mesh.cells, mesh.cell_data["cad:kind"])
                if k[0] == 1
            ]
        )
    )
    r = np.hypot(mesh.points[edge_pts, 0] - 2, mesh.points[edge_pts, 1] - 2)
    assert np.abs(r[r < 1.5] - 1.0).max() < 1e-14
    full = mp.cad_json.read(
        MESHES / "plate_with_hole.cad.json", subdivisions=8, trim=False
    )
    assert sum(len(c.data) for c in full.cells) > sum(len(c.data) for c in mesh.cells)


def test_two_patch_beam_side_channel():
    brep = mp.cad_json.read(MESHES / "two_patch_beam.cad.json").info["brep"]
    assert [f["id"] for f in brep["faces"]] == [2, -1]
    assert brep["faces"][1]["name"] == "right"
    assert brep["faces"][1]["control_point_ids"].size == 9
    assert brep["faces"][1]["loops"][0]["trims"][0]["curve"]["ids"].size == 0
    edges = {e["id"]: e for e in brep["edges"]}
    coupling = edges[10]["topology"]
    assert [
        (r["id"], r["name"], r["trim_index"], r["relative_direction"]) for r in coupling
    ] == [
        (2, "", 1, 1),
        (-1, "right", 3, 0),
    ]
    assert edges[11]["curve"] is None
    assert brep["vertices"][0]["topology"][0]["local_coordinates"] == [2.5, 0.5]
    assert brep["faces"][0]["embedded_edges"][0]["trim_index"] == 4
    assert json.loads(brep["faces"][0]["extra"]) == {
        "embedded_loops": [],
        "embedded_points": [],
    }
    assert json.loads(brep["breps"][1]["extra"]) == {"layer": "beam"}
    assert brep["breps"][1]["name"] == "beam_right"


def test_import_dialect():
    mesh = mp.read(MESHES / "import_square.json", file_format="cad_json")
    brep = mesh.info["brep"]
    # ImportNurbsSbmModeler reads no surface weights: all zero is polynomial.
    assert not brep["faces"][0]["rational"]
    assert all(e["curve"]["rational"] for e in brep["edges"])
    assert json.loads(brep["edges"][4]["extra"]) == {"Layer": "hole"}
    hole = [i for i, c in enumerate(mesh.cells) if c.type == "line3"][0]
    p = mesh.points[np.unique(mesh.cells[hole].data)]
    assert np.abs(np.hypot(p[:, 0] - 0.5, p[:, 1] - 0.5) - 0.25).max() < 1e-15


KRATOS_STYLE = """{
    // Kratos writes control points with ids; files by hand often have none.
    "breps": [{
        "brep_id": 1,
        "faces": [{
            "brep_id": 2,
            "surface": {
                "degrees": [1, 1],
                "knot_vectors": [[0, 1], [0, 1]],
                "control_points": [[[0, 0, 0, 1]], [[1, 0, 0, 1]], [[0, 1, 0, 1]], [[1, 1, 0, 1]]]
            }
        }]
    }]
}"""


@pytest.mark.parametrize("engine", ["core", "python"])
def test_reads_reduced_knots_comments_and_bare_control_points(tmp_path, engine):
    if engine == "core" and _core is None:
        pytest.skip("needs the compiled core")
    path = tmp_path / "hand.cad.json"
    path.write_text(KRATOS_STYLE)
    mesh = (mp.cad_json if engine == "core" else ref).read(path)
    assert [c.type for c in mesh.cells] == ["quad"]
    face = mesh.info["brep"]["faces"][0]
    assert (
        face["reduced_knots"]
        and not face["rational"]
        and face["control_point_ids"].size == 0
    )
    assert mesh.info["brep"]["patches"][0]["knots"][0].tolist() == [0, 0, 1, 1]


BAD = {
    "not_object": "[1]",
    "no_breps": '{"x": 1}',
    "control_mismatch": '{"breps":[{"faces":[{"surface":{"degrees":[1,1],"knot_vectors":[[0,0,1,1],[0,0,1,1]],'
    '"control_points":[[[0,0,0,1]]]}}]}]}',
    "decreasing_knots": '{"breps":[{"faces":[{"surface":{"degrees":[1,1],"knot_vectors":[[1,0],[0,1]],'
    '"control_points":[[[0,0,0,1]],[[1,0,0,1]],[[0,1,0,1]],[[1,1,0,1]]]}}]}]}',
    "zero_weight": '{"breps":[{"faces":[{"surface":{"is_rational":true,"degrees":[1,1],"knot_vectors":[[0,1],[0,1]],'
    '"control_points":[[[0,0,0,0]],[[1,0,0,1]],[[0,1,0,1]],[[1,1,0,1]]]}}]}]}',
    "short_point": '{"breps":[{"edges":[{"3d_curve":{"degree":1,"knot_vector":[0,1],'
    '"control_points":[[1,[0,0,0]],[2,[1,0,0,1]]]}}]}]}',
    "degree_zero": '{"breps":[{"edges":[{"3d_curve":{"degree":0,"knot_vector":[0,1],'
    '"control_points":[[1,[0,0,0,1]]]}}]}]}',
    "mixed_ids": '{"breps":[{"edges":[{"3d_curve":{"degree":1,"knot_vector":[0,1],'
    '"control_points":[[1,[0,0,0,1]],[[1,0,0,1]]]}}]}]}',
    "truncated": '{"breps":[{"faces":[',
    # JSON is UTF-8: a name the core could not hand to Python is refused by
    # both engines (found by mutation fuzzing; the core used to accept it).
    "invalid_utf8": b'{"breps":[{"brep_name":"\xfd"}]}',
    "overlong_utf8": b'{"breps":[{"brep_name":"\xc0\xaf"}]}',
    "surrogate_utf8": b'{"breps":[{"brep_name":"\xed\xa0\x80"}]}',
}


@pytest.mark.parametrize("engine", ["core", "python"])
@pytest.mark.parametrize("name", sorted(BAD))
def test_malformed_files_are_read_errors(tmp_path, engine, name):
    if engine == "core" and _core is None:
        pytest.skip("needs the compiled core")
    path = tmp_path / f"{name}.cad.json"
    body = BAD[name]
    path.write_bytes(body if isinstance(body, bytes) else body.encode())
    reader = _core.cad_json_read if engine == "core" else ref.read
    with pytest.raises(ReadError):
        reader(str(path))


@pytest.mark.parametrize("bad", [-1, 65])
def test_bad_subdivisions_are_value_errors(bad):
    with pytest.raises(ValueError):
        mp.cad_json.read(MESHES / "cylindrical_roof.cad.json", subdivisions=bad)
    with pytest.raises(ValueError):
        ref.read(MESHES / "cylindrical_roof.cad.json", subdivisions=bad)


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", CANONICAL)
def test_rewrite_is_byte_identical_in_both_engines(tmp_path, name):
    mesh = mp.cad_json.read(MESHES / name)
    original = (MESHES / name).read_bytes()
    mp.cad_json.write(tmp_path / "core.cad.json", mesh)
    ref.write(tmp_path / "py.cad.json", mesh)
    assert (tmp_path / "py.cad.json").read_bytes() == original
    assert (tmp_path / "core.cad.json").read_bytes() == original


@needs_core
def test_engines_spell_every_number_alike(tmp_path):
    mesh = mp.cad_json.read(MESHES / "cylindrical_roof.cad.json")
    brep = mesh.info["brep"]
    rng = np.random.default_rng(7)
    cps = rng.standard_normal((6, 3)) * 10.0 ** rng.integers(-30, 30, (6, 3))
    brep["patches"][0]["control_points"] = cps
    brep["patches"][0]["weights"] = rng.uniform(0.1, 3.0, 6)
    mp.cad_json.write(tmp_path / "core.cad.json", mesh, brep=brep)
    ref.write(tmp_path / "py.cad.json", mesh, brep=brep)
    assert (tmp_path / "core.cad.json").read_bytes() == (
        tmp_path / "py.cad.json"
    ).read_bytes()
    back = mp.cad_json.read(tmp_path / "core.cad.json").info["brep"]
    assert np.array_equal(back["patches"][0]["control_points"], cps)


def test_import_dialect_is_written_as_breps(tmp_path):
    mesh = mp.read(MESHES / "import_square.json", file_format="cad_json")
    mp.write(tmp_path / "sq.cad.json", mesh)
    back = mp.read(tmp_path / "sq.cad.json")
    _same_mesh(mesh, back)
    assert "breps" in json.loads((tmp_path / "sq.cad.json").read_text())


def test_plain_mesh_writes_quads_as_bilinear_faces(tmp_path):
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0], [2, 1, 0.5]], float
    )
    mesh = mp.Mesh(
        pts,
        [
            ("quad", [[0, 1, 2, 3], [1, 4, 5, 2]]),
            ("line", [[0, 1]]),
            ("triangle", [[0, 1, 2]]),
        ],
    )
    ref.write(tmp_path / "py.cad.json", mesh)  # warns that the triangle is dropped
    mp.cad_json.write(tmp_path / "core.cad.json", mesh)
    assert (tmp_path / "py.cad.json").read_bytes() == (
        tmp_path / "core.cad.json"
    ).read_bytes()
    back = mp.read(tmp_path / "core.cad.json")
    quads = back.cells[0].data
    assert back.cells[0].type == "quad" and len(quads) == 2
    # Each bilinear face samples back to its quad.
    assert np.allclose(back.points[quads].reshape(-1, 3), pts[[0, 1, 2, 3, 1, 4, 5, 2]])
    with pytest.raises(WriteError):
        mp.cad_json.write(
            tmp_path / "x.cad.json", mp.Mesh(pts, [("triangle", [[0, 1, 2]])])
        )


@pytest.mark.parametrize("engine", ["core", "python"])
def test_inconsistent_brep_is_a_write_error(tmp_path, engine):
    if engine == "core" and _core is None:
        pytest.skip("needs the compiled core")
    mesh = mp.cad_json.read(MESHES / "cylindrical_roof.cad.json")
    brep = mesh.info["brep"]
    brep["faces"][0]["patch"] = 5
    writer = mp.cad_json.write if engine == "core" else ref.write
    with pytest.raises(WriteError):
        writer(tmp_path / "x.cad.json", mesh, brep=brep)


def test_sample_matches_read():
    mesh = mp.cad_json.read(MESHES / "two_patch_beam.cad.json", subdivisions=2)
    _same_mesh(mp.cad_json.sample(mesh.info["brep"], subdivisions=2), mesh)
    _same_mesh(ref.sample(mesh.info["brep"], subdivisions=2), mesh)


def test_strict_core_reads_every_fixture(monkeypatch):
    from meshioplusplus import _fallback

    if _core is None:
        pytest.skip("needs the compiled core")
    _fallback.set_strict_core(True)
    try:
        for name in FIXTURES:
            mp.cad_json.read(MESHES / name)
    finally:
        _fallback.set_strict_core(None)


# ---------------------------------------------------------------------------
# The MDPA sidecar
# ---------------------------------------------------------------------------

_DECK = """Begin ModelPartData
End ModelPartData

Begin Properties 0
End Properties

Begin Nodes
    1 0.0 0.0 0.0
    2 1.0 0.0 0.0
    3 1.0 1.0 0.0
End Nodes

Begin Elements Element2D3N
    1 0 1 2 3
End Elements
"""


def _project(tmp_path, layout):
    deck = tmp_path / "model.mdpa"
    deck.write_text(_DECK)
    geometry = (MESHES / "plate_with_hole.cad.json").read_text()
    if layout == "stem":
        (tmp_path / "model.cad.json").write_text(geometry)
    elif layout == "project_parameters":
        (tmp_path / "geometry.cad.json").write_text(geometry)
        (tmp_path / "ProjectParameters.json").write_text(
            json.dumps(
                {
                    "modelers": [
                        {"modeler_name": "SomethingElse", "Parameters": {}},
                        {
                            "modeler_name": "CadIoModeler",
                            "Parameters": {"geometry_file_name": "geometry"},
                        },
                    ]
                }
            )
        )
    return deck


@pytest.mark.parametrize("layout", ["stem", "project_parameters", "none"])
def test_mdpa_finds_its_cad_json(tmp_path, layout):
    deck = _project(tmp_path, layout)
    found = mp.mdpa.cad_json_path(deck)
    expected = {
        "stem": str(tmp_path / "model.cad.json"),
        "project_parameters": str(tmp_path / "geometry.cad.json"),
        "none": None,
    }[layout]
    assert found == expected
    if _core is not None:
        assert _core.mdpa_cad_json_path(str(deck)) == (expected or "")
    mesh = mp.mdpa.read(deck, cad_json=True)
    assert mesh.cells[0].type == "triangle"
    if layout == "none":
        assert not (mesh.info or {}).get("brep")
    else:
        _same_brep(
            mesh.info["brep"], ref.read_brep(MESHES / "plate_with_hole.cad.json")
        )
    assert not (mp.mdpa.read(deck).info or {}).get("brep")


def test_mdpa_writes_its_cad_json(tmp_path):
    deck = _project(tmp_path, "stem")
    mesh = mp.mdpa.read(deck, cad_json=True)
    out = tmp_path / "out" / "copy.mdpa"
    out.parent.mkdir()
    mp.mdpa.write(out, mesh, cad_json=True)
    assert (out.parent / "copy.cad.json").read_bytes() == (
        MESHES / "plate_with_hole.cad.json"
    ).read_bytes()
    with pytest.raises(ReadError):
        mp.mdpa.read(deck, cad_json=tmp_path / "missing.cad.json")
    with pytest.raises(WriteError):
        mp.mdpa.write(tmp_path / "x.mdpa", mp.mdpa.read(deck), cad_json=True)


# ---------------------------------------------------------------------------
# Kratos itself (optional)
# ---------------------------------------------------------------------------


def _kratos():
    km = pytest.importorskip("KratosMultiphysics")
    pytest.importorskip("KratosMultiphysics.IgaApplication")
    return km


@pytest.mark.parametrize("name", CANONICAL + ["plain"])
def test_kratos_reads_what_meshio_writes(tmp_path, name):
    """Kratos's own CadJsonInput loads the file and evaluates every face where
    meshio++ sampled it (the faces' corner nodes, at their parameters)."""
    km = _kratos()
    if name == "plain":
        pts = np.array(
            [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0], [2, 1, 0.5]], float
        )
        mesh = mp.Mesh(
            pts, [("quad", [[0, 1, 2, 3], [1, 4, 5, 2]]), ("line", [[0, 4]])]
        )
    else:
        mesh = mp.cad_json.read(MESHES / name)
    out = tmp_path / "out.cad.json"
    mp.cad_json.write(out, mesh)
    back = mp.cad_json.read(out, subdivisions=2)
    model = km.Model()
    part = model.CreateModelPart("cad")
    km.CadJsonInput(str(out)).ReadModelPart(part)
    brep = back.info["brep"]
    named = {f["name"]: f for f in brep["faces"] if f["name"]}
    checked = 0
    for f in brep["faces"]:
        geometry = (
            part.GetGeometry(f["id"]) if f["id"] >= 0 else part.GetGeometry(f["name"])
        )
        patch = brep["patches"][f["patch"]]
        (p, q), (ku, kv) = [int(d) for d in patch["degrees"]], patch["knots"]
        for u in np.linspace(ku[p], ku[-p - 1], 5):
            for v in np.linspace(kv[q], kv[-q - 1], 5):
                x = geometry.GlobalCoordinates(km.Vector([float(u), float(v), 0.0]))
                nu, nv = len(ku) - p - 1, len(kv) - q - 1
                ours = ne.surface_point_in_span(
                    p,
                    q,
                    list(ku),
                    list(kv),
                    patch["control_points"].tolist(),
                    patch["weights"].tolist() if f["rational"] else None,
                    nu,
                    ne.find_span(p, list(ku), nu, float(u)),
                    ne.find_span(q, list(kv), nv, float(v)),
                    float(u),
                    float(v),
                )
                assert np.allclose(list(x), ours, rtol=0, atol=1e-12)
                checked += 1
    assert checked and (not named or "right" in named)
