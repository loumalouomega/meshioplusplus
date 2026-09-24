"""Nastran OP2 result files: both engines on pyNastran's test models (real MSC and
NX output, 32- and 64-bit), every step checked against pyNastran's own reading
(frozen in ``pynastran_reference.npz`` by ``tools/gen_nastran_op2_reference.py``)."""

import pathlib
import shutil

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.nastran._nastran import read as read_deck
from meshioplusplus.nastran_op2 import _op2 as py_op2

MESHES = pathlib.Path(__file__).parent / "meshes" / "nastran_op2"
FIXTURES = sorted(MESHES.glob("*.op2"))
# The derived coordinate-system probe has its own reference (cord_reference.npz).
REFERENCE_FIXTURES = [p for p in FIXTURES if not p.stem.endswith("_cord")]
REFERENCE = MESHES / "pynastran_reference.npz"
STATIC = MESHES / "static_solid_shell_bar.op2"


@pytest.fixture(params=["core", "python"])
def engine(request):
    if request.param == "core":
        return meshioplusplus.nastran_op2
    return py_op2


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    for mine, theirs in ((a.point_data, b.point_data), (a.field_data, b.field_data)):
        assert sorted(mine) == sorted(theirs)
        for name in mine:
            np.testing.assert_array_equal(mine[name], theirs[name])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_array_equal(x, y)
    assert sorted(
        (r.kind, r.name, r.dim, r.tag, tuple(r.entries)) for r in a.regions
    ) == sorted((r.kind, r.name, r.dim, r.tag, tuple(r.entries)) for r in b.regions)


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
@pytest.mark.parametrize("step", [0, -1])
def test_engines_agree_on_every_fixture(path, step):
    _same(
        meshioplusplus.nastran_op2.read(path, time_step=step),
        py_op2.read(path, time_step=step),
    )


def _grid_ids(path):
    """GRID id of each point: GEOM1's order, or the sibling deck's."""
    stream, _, tables = py_op2._parse(path.read_bytes())
    ids = py_op2._read_grids(stream, tables)[0]
    if ids:
        return ids
    return read_deck(path.with_suffix(".bdf")).points_id.tolist()


@pytest.mark.parametrize(
    "path", REFERENCE_FIXTURES, ids=[p.name for p in REFERENCE_FIXTURES]
)
def test_every_step_matches_pynastran(engine, path):
    """The done-when: displacements, eigenvectors, SPC/MPC forces, loads,
    temperatures and the centre stress and strain match pyNastran on every step
    (among them a SOL 101 file, ``static_solid_shell_bar``)."""
    ref = np.load(REFERENCE)
    stem = path.stem
    point = {g: i for i, g in enumerate(_grid_ids(path))}
    times = engine.time_values(path)
    rank = {}
    compared = 0
    for step in range(len(times)):
        mesh = engine.read(path, time_step=step)
        sub = int(mesh.field_data["nastran:subcase"][0])
        ana = int(mesh.field_data["nastran:analysis"][0])
        mode = int(mesh.field_data["nastran:mode"][0])
        tag = mode if ana in (2, 8, 9) else rank.setdefault((sub, ana), [0])[0]
        if ana not in (2, 8, 9):
            rank[(sub, ana)][0] += 1
        cell = {}
        for b, eids in enumerate(mesh.cell_data["nastran:eid"]):
            for i, e in enumerate(eids.tolist()):
                cell.setdefault(e, (b, i))
        for key in ref.files:
            parts = key.split("|")
            if parts[-1] in ("ids", "cols") or parts[0] != stem:
                continue
            if "@" in parts[4] or parts[-1] == "gpf":  # multi-valued: tested below
                continue
            ksub, kana, ktag = int(parts[1]), int(parts[2]), int(parts[3])
            if ksub != sub or kana != ana or ktag != tag:
                continue
            name, values, ids = parts[4], ref[key], ref[key + "|ids"]
            if len(parts) == 5:  # nodal
                got = mesh.point_data[name]
                rows = [point[int(g)] for g in ids if int(g) in point]
                sel = [k for k, g in enumerate(ids) if int(g) in point]
                np.testing.assert_allclose(
                    got[rows], values[sel], rtol=1e-6, atol=1e-30
                )
            else:  # element centre values
                got = np.array(
                    [
                        mesh.cell_data[name][cell[int(e)][0]][cell[int(e)][1]]
                        for e in ids
                    ]
                )
                np.testing.assert_allclose(got, values, rtol=1e-6, atol=1e-30)
            compared += len(values)
    if stem not in ("sol401_tstep1",) and times:
        assert compared > 0


def test_sol101_static_reads_the_mesh_and_results(engine):
    mesh = engine.read(STATIC)
    assert {c.type for c in mesh.cells} >= {
        "hexahedron",
        "wedge",
        "tetra",
        "quad",
        "triangle",
        "line",
    }
    assert mesh.field_data["nastran:subcase"][0] == 1
    assert mesh.field_data["nastran:analysis"][0] == 1
    assert mesh.field_data["meshio:time"][0] == 0.0
    for name in (
        "DISPLACEMENT",
        "DISPLACEMENT_ROT",
        "SPC_FORCE",
        "MPC_FORCE",
        "APPLIED_LOAD",
    ):
        assert mesh.point_data[name].shape == (len(mesh.points), 3)
    for name in (
        "STRESS:X",
        "STRESS:VON_MISES",
        "STRESS:X1",
        "STRESS:VON_MISES1",
        "STRAIN:X",
    ):
        assert name in mesh.cell_data
    names = {r.name for r in mesh.regions}
    assert {"PSHELL_4", "PSOLID_2"} <= names


def test_modes_are_steps(engine):
    path = MESHES / "mode_solid_shell_bar.op2"
    times = engine.time_values(path)
    assert len(times) == 3 and times == sorted(times)
    last = engine.read(path, time_step=-1)
    assert last.field_data["nastran:mode"][0] == 3
    assert last.field_data["meshio:time"][0] == times[-1]
    assert "EIGENVECTOR" in last.point_data


def test_64_bit_and_nx_files(engine):
    d173 = engine.read(MESHES / "d173.op2")
    assert len(d173.points) == 172 and "DISPLACEMENT" in d173.point_data
    nx = engine.read(MESHES / "sol401_tstep1.op2", time_step=-1)
    assert len(nx.points) == 9 and {c.type for c in nx.cells} == {"quad"}


def test_thermal_and_sort2(engine):
    mesh = engine.read(MESHES / "time_thermal_elements.op2", time_step=-1)
    assert mesh.point_data["TEMPERATURE"].shape == (len(mesh.points),)
    sort2 = MESHES / "time_thermal_elements_sort2_nx.op2"
    # the SORT2 tables are skipped; the SORT1 twins in the same file are read
    assert len(engine.time_values(sort2)) == 9


def test_without_geometry_the_sibling_deck_gives_the_mesh(engine, tmp_path):
    path = MESHES / "solid_bending_no_geom.op2"
    mesh = engine.read(path)
    full = engine.read(MESHES / "solid_bending.op2")
    assert len(mesh.points) == len(full.points)
    assert sum(len(c.data) for c in mesh.cells) == sum(len(c.data) for c in full.cells)
    by_id = dict(
        zip(
            np.concatenate(full.cell_data["nastran:eid"]).tolist(),
            np.concatenate(full.cell_data["STRESS:VON_MISES"]).tolist(),
        )
    )
    mine = dict(
        zip(
            np.concatenate(mesh.cell_data["nastran:eid"]).tolist(),
            np.concatenate(mesh.cell_data["STRESS:VON_MISES"]).tolist(),
        )
    )
    assert mine == by_id
    alone = tmp_path / "alone.op2"
    shutil.copy(path, alone)
    with pytest.raises(meshioplusplus.ReadError, match="no GEOM1 GRID records"):
        engine.read(alone)


def test_selective_read(engine):
    mesh = engine.read(STATIC, arrays=["DISPLACEMENT", "STRESS:X"])
    assert set(mesh.point_data) == {"DISPLACEMENT"}
    assert "STRESS:X" in mesh.cell_data and "STRESS:Y" not in mesh.cell_data
    bare = engine.read(STATIC, points_only=True)
    assert not bare.point_data and "meshio:time" in bare.field_data


def test_sequence_and_metadata():
    path = MESHES / "mode_solid_shell_bar.op2"
    times = meshioplusplus.read_metadata(path)["time_values"]
    assert _core.read_metadata(str(path))["time_values"] == times
    assert [e["step"] for e in meshioplusplus.sequence_entries(str(path))] == [0, 1, 2]


def test_sniffed_without_the_extension(tmp_path):
    for name in ("static_solid_shell_bar.op2", "d173.op2"):
        path = tmp_path / "results.bin"
        shutil.copy(MESHES / name, path)
        assert meshioplusplus.sniff_format(path) == "nastran_op2"
        assert _core.sniff_format(str(path)) == "nastran_op2"


def test_truncated_file_is_refused(engine, tmp_path):
    path = tmp_path / "cut.op2"
    path.write_bytes(STATIC.read_bytes()[:20000])
    with pytest.raises(meshioplusplus.ReadError):
        engine.read(path)


def test_coordinate_systems_match_pynastran(engine):
    """``static_solid_shell_bar_cord.op2`` puts every GRID of the file in its own
    CORD2R/C/S systems (tools/gen_nastran_cord_reference.py)."""
    ref = np.load(MESHES / "cord_reference.npz")
    mesh = engine.read(MESHES / "static_solid_shell_bar_cord.op2")
    np.testing.assert_allclose(mesh.points, ref["xyz_basic"], rtol=0, atol=1e-12)
    assert np.count_nonzero(mesh.point_data["nastran:cp"]) > 0
    rows = np.searchsorted(ref["nids"], ref["displacement_nids"])
    disp = np.hstack(
        [mesh.point_data["DISPLACEMENT"], mesh.point_data["DISPLACEMENT_ROT"]]
    )
    scale = np.abs(ref["displacement_basic"]).max()
    # the file stores single precision
    np.testing.assert_allclose(
        disp[rows], ref["displacement_basic"], rtol=0, atol=1e-6 * scale
    )


def _cell_lookup(mesh):
    cell = {}
    for b, eids in enumerate(mesh.cell_data["nastran:eid"]):
        for i, e in enumerate(eids.tolist()):
            cell.setdefault(e, (b, i))
    return cell


@pytest.mark.parametrize(
    "path", REFERENCE_FIXTURES, ids=[p.name for p in REFERENCE_FIXTURES]
)
def test_multi_valued_results_match_pynastran(engine, path):
    """Composite plies, plate and solid corners, CBEAM stations and grid point
    forces against pyNastran 1.4.1 (``|cols``: ply id, GRID id, station distance)."""
    ref = np.load(REFERENCE)
    stem = path.stem
    grids = _grid_ids(path)
    point = {g: i for i, g in enumerate(grids)}
    times = engine.time_values(path)
    rank = {}
    compared = 0
    for step in range(len(times)):
        mesh = engine.read(path, time_step=step)
        sub = int(mesh.field_data["nastran:subcase"][0])
        ana = int(mesh.field_data["nastran:analysis"][0])
        mode = int(mesh.field_data["nastran:mode"][0])
        tag = mode if ana in (2, 8, 9) else rank.setdefault((sub, ana), [0])[0]
        if ana not in (2, 8, 9):
            rank[(sub, ana)][0] += 1
        cell = _cell_lookup(mesh)
        for key in ref.files:
            parts = key.split("|")
            if len(parts) != 6 or parts[0] != stem or "@" not in parts[4] + parts[5]:
                if not (len(parts) == 6 and parts[0] == stem and parts[5] == "gpf"):
                    continue
            if (int(parts[1]), int(parts[2]), int(parts[3])) != (sub, ana, tag):
                continue
            name, values = parts[4], ref[key]
            ids, cols = ref[key + "|ids"], ref[key + "|cols"]
            if name.startswith("GRID_FORCE:") and name.count(":") == 2:  # by GRID
                keep = [
                    k for k, g in enumerate(ids.tolist()) if g in point
                ]  # not SPOINTs
                got = mesh.point_data[name][[point[int(ids[k])] for k in keep]]
                np.testing.assert_allclose(got, values[keep], rtol=1e-6, atol=1e-9)
                compared += len(keep)
                continue
            data = mesh.cell_data[name]
            station_sd = mesh.cell_data.get(name.split(":")[0] + ":SD@station")
            for e, col, v in zip(ids.tolist(), cols.tolist(), values.tolist()):
                if int(e) not in cell:  # springs and dampers have no cell
                    continue
                b, i = cell[int(e)]
                if name.endswith("@ply"):
                    k = int(col) - 1
                elif name.endswith("@station"):
                    k = int(np.flatnonzero(station_sd[b][i] == col)[0])
                else:  # a corner or element node: its GRID's place in the cell
                    nodes = mesh.cells[b].data[i].tolist()
                    if (
                        point.get(int(col)) not in nodes
                    ):  # a mid-side node read as linear
                        continue
                    k = nodes.index(point[int(col)])
                np.testing.assert_allclose(data[b][i, k], v, rtol=1e-6, atol=1e-9)
                compared += 1
    if stem in ("static_elements", "static_solid_shell_bar"):
        assert compared > 0
