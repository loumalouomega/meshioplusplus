"""LS-DYNA d3plot state databases: both engines, lasso-python's real test
families, two trimmed LS-DYNA families from Ansys' example data (SPH and
erosion) and the families ``tools/gen_d3plot_reference.py`` writes with
lasso-python, every state checked against lasso-python's own reading (frozen
in ``lasso_reference.npz``); airbags and rigid roads on a file built here."""

import pathlib
import shutil
import struct

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.lsdyna_d3plot import _d3plot as py_d3plot

MESHES = pathlib.Path(__file__).parent / "meshes" / "lsdyna_d3plot"
FAMILIES = sorted(MESHES.glob("*/*/d3plot"))
REFERENCE = MESHES / "lasso_reference.npz"
SHELL_SOLID = MESHES / "generated" / "shell_solid" / "d3plot"


def _id(path):
    return f"{path.parent.parent.name}/{path.parent.name}"


@pytest.fixture(params=["core", "python"])
def engine(request):
    if request.param == "core":
        return meshioplusplus.lsdyna_d3plot
    return py_d3plot


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    for mine, theirs in ((a.point_data, b.point_data), (a.field_data, b.field_data)):
        assert sorted(mine) == sorted(theirs)
        for name in mine:
            assert np.asarray(mine[name]).dtype == np.asarray(theirs[name]).dtype
            np.testing.assert_array_equal(mine[name], theirs[name])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            assert x.dtype == y.dtype
            np.testing.assert_array_equal(x, y)
    assert sorted(
        (r.kind, r.name, r.dim, r.tag, tuple(r.entries)) for r in a.regions
    ) == sorted((r.kind, r.name, r.dim, r.tag, tuple(r.entries)) for r in b.regions)


@pytest.mark.parametrize("path", FAMILIES, ids=[_id(p) for p in FAMILIES])
@pytest.mark.parametrize("step", [0, 1, -1])
def test_engines_agree_on_every_family(path, step):
    n = len(py_d3plot.time_values(path))
    if step >= n:
        pytest.skip("fewer states")
    _same(
        meshioplusplus.lsdyna_d3plot.read(path, time_step=step),
        py_d3plot.read(path, time_step=step),
    )


@pytest.mark.parametrize("path", FAMILIES, ids=[_id(p) for p in FAMILIES])
def test_every_state_matches_lasso(engine, path):
    """The done-when: displacement, plastic strain, stress and deletion match
    lasso-python on every state (the generated shell + solid families delete
    shells and a solid along the way)."""
    ref = np.load(REFERENCE)
    key = _id(path)
    times = engine.time_values(path)
    np.testing.assert_allclose(times, ref[f"{key}:times"], rtol=1e-7)
    coords = ref[f"{key}:coords"].astype(np.float64)
    for step in range(len(times)):
        mesh = engine.read(path, time_step=step)
        np.testing.assert_array_equal(mesh.points, coords)
        if f"{key}:disp" in ref:
            np.testing.assert_array_equal(
                mesh.point_data["displacement"],
                ref[f"{key}:disp"][step].astype(np.float64) - coords,
            )
        for family, cell_types in (
            (
                "solid",
                {
                    "hexahedron",
                    "wedge",
                    "tetra",
                    "pyramid",
                    "tetra10",
                    "hexahedron20",
                    "hexahedron27",
                },
            ),
            ("shell", {"quad", "triangle", "quad8"}),
        ):
            if f"{key}:{family}_ids" not in ref or len(ref[f"{key}:{family}_ids"]) == 0:
                continue
            by_eid = {}
            for b, cells in enumerate(mesh.cells):
                if cells.type not in cell_types:
                    continue
                for i, e in enumerate(mesh.cell_data["lsdyna:eid"][b].tolist()):
                    by_eid[e] = (b, i)
            order = [by_eid[int(e)] for e in ref[f"{key}:{family}_ids"]]
            if f"{key}:{family}_eps" in ref:
                expected = ref[f"{key}:{family}_eps"][step]
                expected = expected.reshape(len(order), -1)
                got = np.array(
                    [mesh.cell_data["effective_plastic_strain"][b][i] for b, i in order]
                ).reshape(len(order), -1)
                np.testing.assert_array_equal(got[:, : expected.shape[1]], expected)
            if f"{key}:{family}_stress" in ref:
                expected = ref[f"{key}:{family}_stress"][step]
                n, p, w = expected.shape
                layout = mesh.field_data.get(
                    "lsdyna_d3plot:layout:stress", np.array([1, 6])
                )
                got = np.array([mesh.cell_data["stress"][b][i] for b, i in order])
                got = got.reshape(n, layout[0], layout[1])[:, :p, :w]
                np.testing.assert_array_equal(got, expected)
            if f"{key}:{family}_alive" in ref:
                expected = (ref[f"{key}:{family}_alive"][step] != 0).astype(np.int8)
                got = np.array([mesh.cell_data["lsdyna:alive"][b][i] for b, i in order])
                np.testing.assert_array_equal(got, expected)


def test_shell_solid_family_deletes_and_yields(engine):
    first = engine.read(SHELL_SOLID, time_step=0)
    last = engine.read(SHELL_SOLID, time_step=-1)
    assert [c.type for c in first.cells] == [
        "hexahedron",
        "wedge",
        "tetra",
        "line",
        "quad",
        "triangle",
    ]
    alive = np.concatenate(last.cell_data["lsdyna:alive"])
    assert alive.dtype == np.int8
    assert alive.sum() == len(alive) - 4  # three shells and the tetrahedron
    eps = np.concatenate([a[:, 0] for a in last.cell_data["effective_plastic_strain"]])
    assert np.nanmax(eps) > 0.0
    assert not np.any(np.concatenate(first.cell_data["lsdyna:alive"]) == 0)
    names = {r.name: r.tag for r in first.regions}
    assert names == {"block": 100, "plate": 200, "struts": 300}
    np.testing.assert_array_equal(
        last.field_data["lsdyna_d3plot:layout:stress"], [3, 6]
    )
    assert last.field_data["meshio:time"][0] == pytest.approx(3.0e-3)


def test_precision_and_splitting_do_not_change_the_reading(engine):
    single = engine.read(SHELL_SOLID, time_step=2)
    split = engine.read(
        MESHES / "generated" / "shell_solid_split" / "d3plot", time_step=2
    )
    double = engine.read(
        MESHES / "generated" / "shell_solid_double" / "d3plot", time_step=2
    )
    _same(single, split)
    # double precision carries more digits of the same (float32-born) values
    np.testing.assert_allclose(
        double.point_data["displacement"],
        single.point_data["displacement"],
        rtol=1e-6,
        atol=1e-7,
    )


def test_sequence_metadata_and_family_glob(tmp_path):
    folder = MESHES / "generated" / "shell_solid_split"
    times = meshioplusplus.read_metadata(folder / "d3plot")["time_values"]
    assert len(times) == 4
    assert _core.read_metadata(str(folder / "d3plot"))["time_values"] == times
    entries = meshioplusplus.sequence_entries(str(folder / "d3plot*"))
    assert [e["step"] for e in entries] == [0, 1, 2, 3]
    assert {pathlib.Path(e["path"]).name for e in entries} == {"d3plot"}


def test_selective_read(engine):
    mesh = engine.read(SHELL_SOLID, arrays=["displacement"], time_step=1)
    assert "displacement" in mesh.point_data and "stress" not in mesh.cell_data
    bare = engine.read(SHELL_SOLID, points_only=True)
    assert "displacement" not in bare.point_data
    assert "meshio:time" in bare.field_data


def test_found_by_name_and_by_content(tmp_path):
    assert meshioplusplus.read(SHELL_SOLID).cells
    path = tmp_path / "crash.bin"
    shutil.copy(SHELL_SOLID, path)
    assert meshioplusplus.sniff_format(path) == "lsdyna_d3plot"
    assert _core.sniff_format(str(path)) == "lsdyna_d3plot"
    big = MESHES / "generated" / "shell_solid_double" / "d3plot"
    assert _core.sniff_format(str(big)) == "lsdyna_d3plot"


def test_a_family_member_is_refused_by_name(engine):
    member = MESHES / "generated" / "shell_solid_split" / "d3plot02"
    with pytest.raises(meshioplusplus.ReadError, match="open the base file"):
        engine.read(member)
    with pytest.raises(meshioplusplus.ReadError, match="open the base file"):
        meshioplusplus.read(member)


def test_truncated_family_is_refused(engine, tmp_path):
    data = SHELL_SOLID.read_bytes()
    path = tmp_path / "d3plot"
    path.write_bytes(data[:300])
    with pytest.raises(meshioplusplus.ReadError):
        engine.read(path)


def test_refused_features_are_named(engine, tmp_path):
    data = bytearray(SHELL_SOLID.read_bytes())
    path = tmp_path / "d3plot"
    femzip = bytearray(data)
    femzip[51 * 4 : 52 * 4] = (76_893_465).to_bytes(4, "little")
    path.write_bytes(bytes(femzip))
    with pytest.raises(meshioplusplus.ReadError, match="femzip"):
        engine.read(path)
    intfor = bytearray(data)
    intfor[11 * 4 : 12 * 4] = (4).to_bytes(4, "little")
    path.write_bytes(bytes(intfor))
    with pytest.raises(meshioplusplus.ReadError, match="intfor"):
        engine.read(path)


def test_d3part_reads_as_a_d3plot(engine, tmp_path):
    """A d3part (file type 5, the parts a *DATABASE_BINARY_D3PART selects)
    has the d3plot's layout."""
    for member in SHELL_SOLID.parent.iterdir():
        data = bytearray(member.read_bytes())
        if member.name == "d3plot":
            data[11 * 4 : 12 * 4] = (5).to_bytes(4, "little")
        (tmp_path / member.name.replace("d3plot", "d3part")).write_bytes(bytes(data))
    part = tmp_path / "d3part"
    assert meshioplusplus.sniff_format(part) == "lsdyna_d3plot"
    for step in (0, -1):
        _same(
            engine.read(part, time_step=step), engine.read(SHELL_SOLID, time_step=step)
        )


BIRD = MESHES / "dyna" / "bird_strike" / "d3plot"
PROJECTILE = MESHES / "dyna" / "projectile" / "d3plot"
QUADRATIC = MESHES / "generated" / "quadratic_rigid" / "d3plot"


def test_sph_particles_match_lasso(engine):
    """SPH particles (a bird strike: 701 particles on composite shells) are
    vertices with their variables, as lasso-python reads them; the element
    deletion flags come before them in a state."""
    ref = np.load(REFERENCE)
    mesh = engine.read(BIRD)
    b = [c.type for c in mesh.cells].index("vertex")
    np.testing.assert_array_equal(
        mesh.cells[b].data[:, 0], ref["dyna/bird_strike:sph_nodes"]
    )
    for ours, theirs in (
        ("sph_radius", "sph_radius"),
        ("sph_pressure", "sph_pressure"),
        ("stress", "sph_stress"),
        ("effective_plastic_strain", "sph_eps"),
        ("density", "sph_density"),
        ("internal_energy", "sph_internal_energy"),
        ("sph_neighbors", "sph_neighbors"),
        ("strain", "sph_strain"),
        ("strain_rate", "sph_strainrate"),
        ("mass", "sph_mass"),
    ):
        expected = ref["dyna/bird_strike:" + theirs][0]
        width = int(np.prod(expected.shape[1:]))
        got = mesh.cell_data[ours][b].reshape(len(expected), -1)[:, :width]
        np.testing.assert_array_equal(got.reshape(expected.shape), expected)
    np.testing.assert_array_equal(
        mesh.cell_data["lsdyna:alive"][b] == 0, ref["dyna/bird_strike:sph_deletion"][0]
    )
    shells = np.concatenate(
        [
            a
            for a, c in zip(mesh.cell_data["lsdyna:alive"], mesh.cells)
            if c.type != "vertex"
        ]
    )
    assert (shells == 1).all()


def test_erosion_of_a_real_run(engine):
    """The projectile family: elements eroded by LS-DYNA (18, then 614)."""
    dead = [
        int(
            (
                np.concatenate(
                    engine.read(PROJECTILE, time_step=k).cell_data["lsdyna:alive"]
                )
                == 0
            ).sum()
        )
        for k in range(2)
    ]
    assert dead == [18, 614]


def test_quadratic_solids_and_rigid_bodies(engine):
    """20- and 27-node hexahedra in VTK's node order (LS-DYNA's), and the
    motion of a rigid body, as lasso-python wrote them."""
    ref = np.load(REFERENCE)
    mesh = engine.read(QUADRATIC, time_step=-1)
    assert [(c.type, c.data.tolist()) for c in mesh.cells] == [
        ("hexahedron20", [list(range(20))]),
        ("hexahedron27", [list(range(20, 47))]),
    ]
    # the edge nodes are the edges' midpoints, the face nodes the faces' centres
    pts = mesh.points
    hex27 = pts[20:47]
    np.testing.assert_allclose(hex27[8], (hex27[0] + hex27[1]) / 2)
    np.testing.assert_allclose(hex27[20], hex27[[0, 3, 7, 4]].mean(0))
    np.testing.assert_allclose(hex27[26], hex27[:8].mean(0))
    assert mesh.field_data["lsdyna:rigid_body_part"].tolist() == [20]
    for ours, theirs in (
        ("coordinates", "rigid_coordinates"),
        ("rotation", "rigid_rotation"),
        ("velocity", "rigid_velocity"),
        ("rotational_velocity", "rigid_rot_velocity"),
        ("acceleration", "rigid_acceleration"),
        ("rotational_acceleration", "rigid_rot_acceleration"),
    ):
        np.testing.assert_array_equal(
            mesh.field_data["lsdyna:rigid_body_" + ours],
            ref["generated/quadratic_rigid:" + theirs][-1],
        )


def _particles_family():
    """One shell on nodes 1-4 (part 1), SPH particles on nodes 5 and 6
    (material 2), an airbag (id 77) of two particles, a one-segment rigid road
    (id 9) and a reduced rigid body: NDIM 9, two states, laid out as the
    database manual (and lasso-python's reader) says. lasso-python cannot
    write airbags or roads, and no real file with them was available."""
    out = []

    def i(v):
        out.append(struct.pack("<i", v))

    def f(v):
        out.append(struct.pack("<f", v))

    head = {11: 1, 15: 9, 16: 6, 17: 6, 20: 1, 31: 1, 32: 1, 33: 7, 36: -10001}
    head.update({37: 2, 43: 1000, 44: 1000, 51: 2, 54: 1})
    for k in range(64):
        if k < 10:
            out.append(b"    ")
        elif k == 14:
            f(971.0)
        else:
            i(head.get(k, 0))
    for v in (11, 1, 1, 6, 1, 1, 1, 1, 0, 1, 0):  # SPH flags
        i(v)
    for v in (4, 5, 2, 2):  # airbag: geometry, particle, bag variables, particles
        i(v)
    names = ["Start", "Npart", "Bag ID", "NGas", "GasC ID", "Pos x", "Pos y"]
    names += ["Pos z", "Mass", "Act Gas", "Bag Vol"]
    for t in (1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 2):
        i(t)
    for name in names:
        for c in name.ljust(8)[:8]:
            i(ord(c))
    coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (5, 0, 0), (6, 0, 0)]
    for c in coords:
        for x in c:
            f(x)
    for v in (1, 2, 3, 4, 1):  # the shell
        i(v)
    for v in (1, 2, 1, 5, 0):  # a rigid body: part 2, node 5, no active node
        i(v)
    for v in (5, 2, 6, 2):  # SPH (node, material)
        i(v)
    for v in (1, 2, 77, 1):  # airbag 77: particles 1..2
        i(v)
    for v in (4, 1, 1, 1, 101, 102, 103, 104):  # road nodes
        i(v)
    for c in ((0, 0, -1), (2, 0, -1), (2, 2, -1), (0, 2, -1)):
        for x in c:
            f(x)
    for v in (9, 1, 101, 102, 103, 104):  # road 9, one segment
        i(v)
    f(-999999.0)
    for k in range(2):
        f(0.5 * k)
        for c in coords:
            f(c[0])
            f(c[1])
            f(c[2] + 0.1 * k)
        for c in range(6):
            f(10.0 * k + c)
        f(0.25 * k)
        f(1.0)  # the shell lives
        for p in range(2):  # SPH: material (negative once deleted), variables
            f(-2.0 if (k == 1 and p == 1) else 2.0)
            for v in range(13):
                f(100.0 * k + 10.0 * p + v)
        i(2 if k == 0 else 1)  # airbag: active gas (an integer), volume
        f(3.0 + k)
        for p in range(2):  # particles: gas id (an integer), x, y, z, mass
            i(1)
            f(10.0 + p)
            f(20.0 + k)
            f(30.0)
            f(0.5 + p)
        for v in (0.0, 0.0, -0.1 * k, 0.0, 0.0, -1.0):  # road: displacement, velocity
            f(v)
        for v in range(12):  # rigid body: coordinates, rotation
            f(k + v / 10.0)
    return b"".join(out)


def test_airbags_rigid_roads_and_sph_deletion(engine, tmp_path):
    path = tmp_path / "d3plot"
    path.write_bytes(_particles_family())
    mesh = engine.read(path, time_step=1)
    assert [(c.type, c.data.tolist()) for c in mesh.cells] == [
        ("quad", [[0, 1, 2, 3]]),
        ("vertex", [[4], [5]]),  # SPH
        ("vertex", [[6], [7]]),  # airbag particles
        ("quad", [[8, 9, 10, 11]]),  # the road
    ]
    np.testing.assert_allclose(mesh.points[6:8], [[10, 21, 30], [11, 21, 30]])
    np.testing.assert_allclose(mesh.points[8], [0, 0, -1])
    assert np.isnan(mesh.point_data["displacement"][6:]).all()
    assert [a.tolist() for a in mesh.cell_data["lsdyna:alive"]] == [
        [1],
        [1, 0],
        [1, 1],
        [1],
    ]
    assert [a.tolist() for a in mesh.cell_data["lsdyna:part"]] == [
        [1],
        [2, 2],
        [77, 77],
        [9],
    ]
    np.testing.assert_array_equal(mesh.cell_data["sph_radius"][1], [100.0, 110.0])
    np.testing.assert_array_equal(mesh.cell_data["airbag_gasc_id"][2], [1.0, 1.0])
    np.testing.assert_array_equal(mesh.cell_data["airbag_mass"][2], [0.5, 1.5])
    assert mesh.field_data["lsdyna:airbag:act_gas"].tolist() == [1.0]
    assert mesh.field_data["lsdyna:airbag:bag_vol"].tolist() == [4.0]
    np.testing.assert_allclose(mesh.field_data["lsdyna:road_velocity"], [[0, 0, -1]])
    np.testing.assert_allclose(
        mesh.field_data["lsdyna:rigid_body_coordinates"], [[1.0, 1.1, 1.2]], rtol=1e-6
    )
    assert "lsdyna:rigid_body_velocity" not in mesh.field_data  # reduced (NDIM 9)
    names = {r.name: r.entries.tolist() for r in mesh.regions}
    assert names == {
        "Part 1": [0],
        "Part 2": [1, 2],
        "Airbag 77": [3, 4],
        "Rigid road 9": [5],
    }
