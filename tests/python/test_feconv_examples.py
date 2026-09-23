"""FEconv's example meshes, read from a local checkout.

FEconv (https://github.com/victorsndvg/FEconv) ships small meshes in many
formats, several of them the same mesh written by different tools. They are
GPL-3 and not part of this repository, so this module is skipped unless
``MESHIOPLUSPLUS_FECONV_DIR`` points at a checkout (or at its ``examples``)::

    git clone --depth 1 https://github.com/victorsndvg/FEconv
    MESHIOPLUSPLUS_FECONV_DIR=$PWD/FEconv pytest tests/python/test_feconv_examples.py

For every file in a format both engines read, it checks that the C++ core and
the Python reader agree, that every solid is positively oriented and that every
mid-edge node lies near its edge; then that the files known to hold the same
mesh read as the same cells. The quirks these files exposed are reproduced by
the committed fixtures of ``tools/gen_feconv_quirk_fixtures.py``.
"""

import importlib
import os
import pathlib
import warnings

import numpy as np
import pytest
from meshioplusplus import _core
from meshioplusplus._exceptions import ReadError


def _root():
    d = os.environ.get("MESHIOPLUSPLUS_FECONV_DIR")
    if not d:
        return None
    d = pathlib.Path(d)
    return d / "examples" if (d / "examples").is_dir() else d


ROOT = _root()
pytestmark = pytest.mark.skipif(
    ROOT is None or not ROOT.is_dir(),
    reason="set MESHIOPLUSPLUS_FECONV_DIR to a FEconv checkout",
)

_CORE = {
    "flux": _core.flux_read,
    "vtu": _core.vtu_read,
    "medit": _core.medit_read_ascii,
    "gmsh": _core.gmsh_read,
    "ansys": _core.ansys_read,
    "mphtxt": _core.mphtxt_read,
    "unv": _core.unv_read,
}

# Files that are not meshes, or whose cells are inverted in the file itself.
_NO_MESH = {"test/mphtxt/geo6.mphtxt"}  # a COMSOL geometry, no Mesh object
_INVERTED_IN_FILE = {
    # Fluent's UNV export: every tetrahedron negative, under the descriptor 111
    # whose Salome and I-DEAS files read positive.
    "test/unv/malla_cfd_fluent_choose_interior.unv",
}


def _format(path):
    suffix = path.suffix.lower()
    if suffix == ".msh":
        head = path.read_bytes()[:256].lstrip()
        if head.startswith(b"$"):
            return "gmsh"
        return "ansys" if head.startswith(b"(") else None
    return {
        ".pf3": "flux",
        ".vtu": "vtu",
        ".mesh": "medit",
        ".mphtxt": "mphtxt",
        ".unv": "unv",
    }.get(suffix)


def _files():
    if ROOT is None or not ROOT.is_dir():
        return []
    return [
        p.relative_to(ROOT).as_posix()
        for p in sorted(ROOT.rglob("*"))
        if p.is_file() and p.stat().st_size > 0 and _format(p)
    ]


def _read(rel, engine):
    path = ROOT / rel
    fmt = _format(path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if engine == "core":
            return _CORE[fmt](str(path))
        return importlib.import_module(f"meshioplusplus.{fmt}")._py_read(path)


_CORNERS = {
    "tetra": (0, 1, 2, 3),
    "tetra10": (0, 1, 2, 3),
    "pyramid": (0, 1, 3, 4),
    "wedge": (0, 1, 2, 3),
    "wedge15": (0, 1, 2, 3),
    "hexahedron": (0, 1, 3, 4),
    "hexahedron20": (0, 1, 3, 4),
    "hexahedron27": (0, 1, 3, 4),
}
_MIDS = {
    "line3": [(2, 0, 1)],
    "triangle6": [(3, 0, 1), (4, 1, 2), (5, 2, 0)],
    "quad8": [(4, 0, 1), (5, 1, 2), (6, 2, 3), (7, 3, 0)],
    "quad9": [(4, 0, 1), (5, 1, 2), (6, 2, 3), (7, 3, 0)],
    "tetra10": [(4, 0, 1), (5, 1, 2), (6, 2, 0), (7, 0, 3), (8, 1, 3), (9, 2, 3)],
}
_HEX_MIDS = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
_HEX_MIDS += [(0, 4), (1, 5), (2, 6), (3, 7)]
_MIDS["hexahedron20"] = [(8 + k, a, b) for k, (a, b) in enumerate(_HEX_MIDS)]
_MIDS["hexahedron27"] = _MIDS["hexahedron20"]


def _points3(mesh):
    p = np.asarray(mesh.points, dtype=float)
    return (
        p
        if p.shape[1] == 3
        else np.column_stack([p, np.zeros((len(p), 3 - p.shape[1]))])
    )


@pytest.mark.parametrize("rel", _files())
def test_engines_agree_and_cells_are_well_formed(rel):
    if rel in _NO_MESH:
        for engine in ("core", "python"):
            with pytest.raises(Exception):
                _read(rel, engine)
        return
    a, b = _read(rel, "core"), _read(rel, "python")
    np.testing.assert_allclose(a.points, b.points)
    assert [(c.type, len(c.data)) for c in a.cells] == [
        (c.type, len(c.data)) for c in b.cells
    ]
    for x, y in zip(a.cells, b.cells):
        if x.type.startswith("polyhedron"):
            assert [[list(f) for f in c] for c in x.data] == [
                [list(f) for f in c] for c in y.data
            ]
        else:
            np.testing.assert_array_equal(x.data, y.data)

    pts = _points3(a)
    for block in a.cells:
        if block.type.startswith("poly") or not len(block.data):
            continue
        data = np.asarray(block.data)
        if block.type in _CORNERS and rel not in _INVERTED_IN_FILE:
            i, j, k, m = _CORNERS[block.type]
            x = pts[data]
            vol = np.einsum(
                "ij,ij->i",
                np.cross(x[:, j] - x[:, i], x[:, k] - x[:, i]),
                x[:, m] - x[:, i],
            )
            assert (vol > 0).all(), f"{block.type}: {(vol <= 0).sum()} inverted"
        for mid, i, j in _MIDS.get(block.type, []):
            # curved meshes move mid-edge nodes, but not by a quarter of the edge
            dev = np.linalg.norm(
                pts[data[:, mid]] - (pts[data[:, i]] + pts[data[:, j]]) / 2, axis=1
            )
            edge = np.linalg.norm(pts[data[:, i]] - pts[data[:, j]], axis=1)
            assert (
                dev <= 0.25 * edge
            ).all(), f"{block.type}: mid-edge node {mid} misplaced"


# ---------------------------------------------------------------------------
# Files holding the same mesh
# ---------------------------------------------------------------------------


def _cell_sets(mesh, ref_points, dims):
    """The cells of the given dimensions as sets of reference-point ids (by
    nearest coordinates: the writers round differently), corners only. The
    meshes are compared up to a translation (one copy of the cube sits a unit
    higher)."""
    pts = _points3(mesh)
    ref = _points3(ref_points)
    pts = pts - pts.min(axis=0) + ref.min(axis=0)
    scale = float(np.ptp(ref, axis=0).max())
    corners = {"tetra": 4, "tetra10": 4, "hexahedron": 8, "hexahedron20": 8}
    corners.update({"hexahedron27": 8, "wedge": 6, "triangle": 3, "triangle6": 3})
    corners.update({"quad": 4, "quad8": 4, "quad9": 4})
    dim = {"tetra": 3, "tetra10": 3, "hexahedron": 3, "hexahedron20": 3}
    dim.update({"hexahedron27": 3, "wedge": 3, "triangle": 2, "triangle6": 2})
    dim.update({"quad": 2, "quad8": 2, "quad9": 2})
    rows = [
        np.asarray(block.data)[:, : corners[block.type]]
        for block in mesh.cells
        if dim.get(block.type) in dims and len(block.data)
    ]
    ids = {}
    for p in np.unique(np.concatenate([r.ravel() for r in rows])) if rows else []:
        d = np.linalg.norm(ref - pts[p], axis=1)
        ids[int(p)] = int(d.argmin())
        assert d.min() < 1e-4 * scale
    return {frozenset(ids[int(v)] for v in c) for r in rows for c in r}


def _twin(a, b, dims=(3,)):
    if not (ROOT / a).is_file() or not (ROOT / b).is_file():
        pytest.skip(f"{a} or {b} not in this checkout")
    ma, mb = _read(a, "core"), _read(b, "core")
    sa = _cell_sets(ma, mb, dims)
    sb = _cell_sets(mb, mb, dims)
    assert sa and sa == sb


# The same 69-point, 224-tetrahedron cube in every format FEconv converts to.
_CUBE = [
    "mesh.mesh",
    "mesh.mphtxt",
    "mesh.msh",
    "ansys_mesh.msh",
    "mesh.pf3",
    "mesh.vtu",
    "mesh_.vtu",
    "test/vtu/fluid 1_1.vtu",
    "test/vtu/solid 1_1.vtu",
]


@pytest.mark.parametrize("rel", _CUBE)
def test_the_cube_is_the_same_everywhere(rel):
    _twin(rel, "mesh.unv")


_MPHTXT_UNV = [
    ("2solidcubes", (3,)),
    ("2squarefaces", (2,)),
    ("4quads", (2,)),
    ("hexacubelimite", (3,)),
    ("hexap2", (3,)),
    ("isogrid-mesh", (3,)),
    ("mesh-geo8", (2,)),
    ("quadp2", (2,)),
    ("squarefecube", (2,)),
    ("surfacesphere", (2,)),
    ("tetrap2", (3,)),
    ("triap2", (2,)),
]


@pytest.mark.parametrize("name, dims", _MPHTXT_UNV)
def test_comsol_files_match_their_unv_twins(name, dims):
    _twin(f"test/mphtxt/{name}.mphtxt", f"test/unv/{name}.unv", dims)


@pytest.mark.parametrize(
    "pf3, unv, dims",
    [
        ("EXAMPLE", "EXAMPLE", (3,)),
        ("MALLA3D", "MALLA3D", (3,)),
        ("RECTANGULO", "RECTANGILO", (2,)),
    ],
)
def test_flux_files_match_their_unv_twins(pf3, unv, dims):
    _twin(f"test/pf3/{pf3}.pf3", f"test/unv/{unv}.unv", dims)


def test_the_two_object_cubes_are_two_cubes():
    # The reviewer's report: two COMSOL Mesh objects, 12 tetrahedra each.
    rel = "test/mphtxt/2objectcubes.mphtxt"
    if not (ROOT / rel).is_file():
        pytest.skip(f"{rel} not in this checkout")
    mesh = _read(rel, "core")
    tets = [np.asarray(c.data) for c in mesh.cells if c.type == "tetra"]
    assert [len(t) for t in tets] == [12, 12]
    pts = _points3(mesh)
    for t in tets:
        x = pts[t]
        vol = np.einsum(
            "ij,ij->i",
            np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]),
            x[:, 3] - x[:, 0],
        )
        assert np.isclose(vol.sum() / 6, 1.0)
    assert sorted(r.name for r in mesh.regions) == ["mesh1", "mesh2"]


def test_the_files_the_old_readers_refused_read():
    for rel in [
        "test/quad-gmsh.msh",
        "gmsh_mesh.msh",
        "mesh.mesh",
        "test/surfacesphere.mphtxt",
        "test/msh/triamesh2.msh",
        "test/unv/salida_code_aster_series.unv",
        "test/unv/unvmesh2414.unv",
        "test/pf3/MIXED-PENTAHEDRON-HEXAHEDRON-PYRAMID_P2.PF3",
    ]:
        if (ROOT / rel).is_file():
            try:
                mesh = _read(rel, "core")
            except ReadError as exc:
                pytest.fail(f"{rel}: {exc}")
            assert len(mesh.points) > 0, rel
