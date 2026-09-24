"""MSC Nastran HDF5 results (``.h5``, read-only): both engines against real MSC output.

The files under ``meshes/nastran_h5`` were written by MSC Nastran 2020 (see the
README there); ``pynastran_reference.npz`` is pyNastran's own reading of them
(``tools/gen_nastran_h5_reference.py``). The edge cases no real file covers
(quadratic solids, CP != 0, other vendors, broken tables) are built here with h5py
in the MSC schema.
"""

import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core, _sequence
from meshioplusplus._exceptions import ReadError, WriteError
from meshioplusplus._fallback import set_strict_core
from meshioplusplus.nastran_h5 import _nastran_h5 as py_nh5

h5py = pytest.importorskip("h5py")

FIXTURES = pathlib.Path(__file__).parent / "meshes" / "nastran_h5"
NAMES = sorted(p.stem for p in FIXTURES.glob("*.h5"))

pytestmark = pytest.mark.skipif(
    not getattr(_core, "__has_hdf5__", False), reason="core built without HDF5"
)

# Edge (a, b) of each mid-node, in meshio++ order, for a straight-sided element.
MIDS = {
    "hexahedron20": [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6)]
    + [(6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)],
    "wedge15": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
    "tetra10": [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
    "pyramid13": [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)],
    "triangle6": [(0, 1), (1, 2), (2, 0)],
    "quad8": [(0, 1), (1, 2), (2, 3), (3, 0)],
}


def path_of(name):
    path = FIXTURES / f"{name}.h5"
    with open(path, "rb") as f:
        assert not f.read(24).startswith(
            b"version https://git-lfs"
        ), f"{path.name} is an unfetched Git-LFS pointer; run `git lfs pull`"
    return str(path)


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read signature."""
    if request.param == "core":
        # the C++ reader itself, not the wrapper that falls back to Python
        def read(path, points_only=False, arrays=None, time_step=0):
            return _core.nastran_h5_read(str(path), points_only, arrays, time_step)

        return read
    return py_nh5.read


def signed_volumes(mesh, block):
    """Sum of the signed tetra volumes of a linear-corner decomposition."""
    dec = {
        "tetra": [(0, 1, 2, 3)],
        "pyramid": [(0, 1, 2, 4), (0, 2, 3, 4)],
        "wedge": [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)],
        "hexahedron": [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7)]
        + [(1, 3, 4, 6)],
    }[block.type.rstrip("0123456789")]
    p = mesh.points
    c = np.asarray(block.data)
    vol = 0.0
    for a, b, d, e in dec:
        vol = (
            vol
            + np.einsum(
                "ij,ij->i",
                np.cross(p[c[:, b]] - p[c[:, a]], p[c[:, d]] - p[c[:, a]]),
                p[c[:, e]] - p[c[:, a]],
            )
            / 6.0
        )
    return vol


# ---------------------------------------------------------------------------
# A writer of small files in the MSC schema, for what no real file covers.
# ---------------------------------------------------------------------------


def _table(group, name, rows, fields):
    """A compound dataset; ``fields`` is [(name, dtype, shape or ())]."""
    dt = np.dtype([(n, t, s) if s else (n, t) for n, t, s in fields])
    data = np.zeros(len(rows), dtype=dt)
    for i, row in enumerate(rows):
        for (n, _, _), v in zip(fields, row):
            data[n][i] = v
    group.create_dataset(name, data=data)


def write_msc(
    path,
    grids,
    elements=(),
    nodal=None,
    domains=None,
    version=b"msc20200",
    cp=None,
    spoints=(),
    index=True,
):
    """grids: {id: xyz}; elements: [(card, width, [(eid, pid, [g...])])];
    nodal: {table: {domain: [(id, x, y, z, rx, ry, rz)]}}; domains: [(id, time, mode)].
    """
    with h5py.File(path, "w") as f:
        f.attrs["SCHEMA"] = np.array([20200])
        nas = f.create_group("NASTRAN")
        if version is not None:
            nas.attrs["VERSION"] = np.bytes_(version)
        node = nas.create_group("INPUT/NODE")
        ids = sorted(grids)
        _table(
            node,
            "GRID",
            [(g, (cp or {}).get(g, 0), grids[g], 0, 0, 0, 1) for g in ids],
            [
                ("ID", "<i8", ()),
                ("CP", "<i8", ()),
                ("X", "<f8", (3,)),
                ("CD", "<i8", ()),
            ]
            + [("PS", "<i8", ()), ("SEID", "<i8", ()), ("DOMAIN_ID", "<i8", ())],
        )
        if spoints:
            _table(
                node,
                "SPOINT",
                [(s, 1) for s in spoints],
                [("ID", "<i8", ())] + [("DOMAIN_ID", "<i8", ())],
            )
        el = nas.create_group("INPUT/ELEMENT")
        for card, width, rows in elements:
            _table(
                el,
                card,
                [(e, p, list(g) + [0] * (width - len(g)), 1) for e, p, g in rows],
                [("EID", "<i8", ()), ("PID", "<i8", ()), ("G", "<i8", (width,))]
                + [("DOMAIN_ID", "<i8", ())],
            )
        if nodal:
            res = nas.create_group("RESULT")
            _table(
                res,
                "DOMAINS",
                [(d, 1, 0, 1, t, 0.0, m) for d, t, m in domains],
                [("ID", "<i8", ()), ("SUBCASE", "<i8", ()), ("STEP", "<i8", ())]
                + [("ANALYSIS", "<i8", ()), ("TIME_FREQ_EIGR", "<f8", ())]
                + [("EIGI", "<f8", ()), ("MODE", "<i8", ())],
            )
            ng = res.create_group("NODAL")
            ig = f.create_group("INDEX/NASTRAN/RESULT/NODAL")
            for table, per_domain in nodal.items():
                rows, idx = [], []
                for d, drows in per_domain.items():
                    idx.append((d, len(rows), len(drows)))
                    rows += [tuple(r) + (d,) for r in drows]
                _table(
                    ng,
                    table,
                    rows,
                    [("ID", "<i8", ())]
                    + [(c, "<f8", ()) for c in ("X", "Y", "Z", "RX", "RY", "RZ")]
                    + [("DOMAIN_ID", "<i8", ())],
                )
                if index:
                    _table(
                        ig,
                        table,
                        idx,
                        [
                            ("DOMAIN_ID", "<i8", ()),
                            ("POSITION", "<i8", ()),
                            ("LENGTH", "<i8", ()),
                        ],
                    )
    return str(path)


UNIT_HEX = {
    1: (0, 0, 0),
    2: (1, 0, 0),
    3: (1, 1, 0),
    4: (0, 1, 0),
    5: (0, 0, 1),
    6: (1, 0, 1),
    7: (1, 1, 1),
    8: (0, 1, 1),
}


def quadratic_model(tmp_path):
    """A unit CHEXA20, CPENTA15, CTETRA10, CPYRAM13, CTRIA6 and CQUAD8 in Nastran order."""
    grids = dict(UNIT_HEX)
    nid = [100]

    def mid(a, b):
        nid[0] += 1
        grids[nid[0]] = tuple((np.add(grids[a], grids[b]) / 2).tolist())
        return nid[0]

    def mids(edges):
        return [mid(a, b) for a, b in edges]

    # Nastran mid-side order: bottom ring, the vertical edges, top ring.
    hexa = [1, 2, 3, 4, 5, 6, 7, 8] + mids(
        [(1, 2), (2, 3), (3, 4), (4, 1), (1, 5), (2, 6), (3, 7), (4, 8)]
        + [(5, 6), (6, 7), (7, 8), (8, 5)]
    )
    penta = [1, 2, 4, 5, 6, 8] + mids(
        [(1, 2), (2, 4), (4, 1), (1, 5), (2, 6), (4, 8), (5, 6), (6, 8), (8, 5)]
    )
    tetra = [1, 2, 4, 5] + mids([(1, 2), (2, 4), (4, 1), (1, 5), (2, 5), (4, 5)])
    grids[9] = (0.5, 0.5, 2.0)
    pyram = [5, 6, 7, 8, 9] + mids(
        [(5, 6), (6, 7), (7, 8), (8, 5), (5, 9), (6, 9), (7, 9), (8, 9)]
    )
    tria = [1, 2, 3] + mids([(1, 2), (2, 3), (3, 1)])
    quad = [1, 2, 3, 4] + mids([(1, 2), (2, 3), (3, 4), (4, 1)])
    return write_msc(
        tmp_path / "quadratic.h5",
        grids,
        [
            ("CHEXA", 20, [(1, 1, hexa)]),
            ("CPENTA", 15, [(2, 1, penta)]),
            ("CPYRAM", 13, [(3, 1, pyram)]),
            ("CQUAD8", 8, [(4, 2, quad)]),
            ("CTETRA", 10, [(5, 1, tetra)]),
            ("CTRIA6", 6, [(6, 2, tria)]),
        ],
    )


def hex_with_results(tmp_path, **kw):
    rows = {
        1: [(g, g, 0, 0, 0, 0, 0) for g in UNIT_HEX],
        2: [(g, -g, 0, 0, 0, 0, 1) for g in UNIT_HEX if g != 3],
    }
    return write_msc(
        tmp_path / "hex.h5",
        UNIT_HEX,
        [("CHEXA", 20, [(7, 3, list(UNIT_HEX))])],
        nodal={"DISPLACEMENT": rows},
        domains=[(1, 0.5, 0), (2, 1.5, 0)],
        **kw,
    )


# ---------------------------------------------------------------------------
# The real files
# ---------------------------------------------------------------------------


class TestRealFiles:
    def test_static_model(self, engine):
        mesh = engine(path_of("static_elements"))
        assert mesh.points.shape == (40, 3)
        types = [c.type for c in mesh.cells]
        assert {
            "hexahedron",
            "wedge",
            "tetra",
            "quad",
            "triangle",
            "line",
            "vertex",
        } == set(types)
        eids = np.concatenate(mesh.cell_data["nastran:eid"])
        assert len(eids) == sum(len(c.data) for c in mesh.cells) == 34
        pids = np.concatenate(mesh.cell_data["nastran:pid"])
        assert (pids[types_of(mesh) == "vertex"] == -1).all()  # CONM2 has no property

    def test_regions_are_named_after_their_property_card(self, engine):
        mesh = engine(path_of("static_elements"))
        names = {r.name: r for r in mesh.regions}
        assert {"PSHELL_4", "PBAR_1", "PROD_3", "PSHEAR_8"} <= set(names)
        eids = np.concatenate(mesh.cell_data["nastran:eid"])
        pids = np.concatenate(mesh.cell_data["nastran:pid"])
        shell = names["PSHELL_4"]
        assert shell.dim == 2 and shell.tag == 4
        assert sorted(eids[shell.entries]) == sorted(eids[pids == 4])

    @pytest.mark.parametrize("name", NAMES)
    def test_solids_are_not_inverted(self, engine, name):
        mesh = engine(path_of(name))
        for block in mesh.cells:
            if block.type in ("tetra", "pyramid", "wedge", "hexahedron"):
                assert (signed_volumes(mesh, block) > 0).all(), block.type

    def test_steps_of_each_solution(self):
        assert _sequence.num_steps(path_of("static_elements")) == 1
        assert _sequence.num_steps(path_of("modes_elements")) == 3
        assert _sequence.num_steps(path_of("freq_elements")) == 5
        # static (domain 1) + 4 buckling eigenvectors; DOMAINS rows no table uses are not steps
        assert _sequence.num_steps(path_of("buckling_solid_shell_bar")) == 9
        assert meshioplusplus.read_metadata(path_of("time_thermal_elements"))[
            "time_values"
        ] == pytest.approx([0.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0])
        assert meshioplusplus.read_metadata(path_of("freq_elements"))[
            "time_values"
        ] == pytest.approx([1e-5, 10.0, 20.0, 30.0, 40.0])

    def test_a_mode_step(self, engine):
        mesh = engine(path_of("modes_elements"), time_step=2)
        fd = {k: v.tolist() for k, v in mesh.field_data.items()}
        assert fd["nastran:mode"] == [3] and fd["nastran:subcase"] == [1]
        assert fd["nastran:analysis"] == [2] and fd["nastran:domain"] == [6]
        assert mesh.point_data["EIGENVECTOR"].shape == (40, 3)
        assert mesh.point_data["EIGENVECTOR_ROT"].shape == (40, 3)
        assert not np.isnan(mesh.point_data["EIGENVECTOR"]).any()

    def test_spoint_rows_are_not_points(self, engine):
        # EIGENVECTOR has 43 rows per mode: the 40 GRIDs and 3 SPOINTs
        with h5py.File(path_of("modes_elements")) as f:
            rows = f["/INDEX/NASTRAN/RESULT/NODAL/EIGENVECTOR"][0]["LENGTH"]
        assert rows == 43
        assert (
            engine(path_of("modes_elements")).point_data["EIGENVECTOR"].shape[0] == 40
        )

    def test_complex_tables_are_real_and_imaginary_pairs(self, engine):
        mesh = engine(path_of("freq_elements"), time_step=-1)
        for base in ("DISPLACEMENT", "DISPLACEMENT_ROT", "SPC_FORCE"):
            assert mesh.point_data[base + "_real"].shape == (40, 3)
            assert mesh.point_data[base + "_imag"].shape == (40, 3)
        assert mesh.field_data["meshio:time"][0] == pytest.approx(40.0)

    def test_scalar_nodal_table(self, engine):
        mesh = engine(path_of("time_thermal_elements"), time_step=-1)
        assert mesh.point_data["TEMPERATURE"].shape == (9,)

    def test_tables_can_cover_different_domains(self, engine):
        path = path_of("buckling_solid_shell_bar")
        static = engine(path, time_step=0)
        mode = engine(path, time_step=1)
        assert "DISPLACEMENT" in static.point_data
        assert "EIGENVECTOR" not in static.point_data
        assert "EIGENVECTOR" in mode.point_data
        assert "DISPLACEMENT" not in mode.point_data
        assert "SPC_FORCE" in static.point_data and "SPC_FORCE" in mode.point_data
        assert mode.field_data["meshio:time"][0] == pytest.approx(-4.93586081e10)

    def test_element_results_take_the_centre_value(self, engine):
        path = path_of("static_elements")
        mesh = engine(path)
        with h5py.File(path) as f:
            hexa = f["/NASTRAN/RESULT/ELEMENTAL/STRESS/HEXA"][0]
            tria = f["/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3"][()]
        eids = np.concatenate(mesh.cell_data["nastran:eid"])
        types = types_of(mesh)
        sx = np.concatenate(mesh.cell_data["STRESS:X"])
        cell = np.flatnonzero((eids == hexa["EID"]) & (types == "hexahedron"))[0]
        assert sx[cell] == hexa["X"][0]  # index 0 is the centre ("CEN")
        x1 = np.concatenate(mesh.cell_data["STRESS:X1"])
        for row in tria:
            cell = np.flatnonzero((eids == row["EID"]) & (types == "triangle"))[0]
            assert x1[cell] == row["X1"]
        # a line element has no STRESS:X
        assert np.isnan(sx[types == "line"]).all()

    def test_per_ply_and_per_grid_tables_are_skipped(self, engine):
        mesh = engine(path_of("static_elements"))
        assert not any(k.startswith("GRID_FORCE") for k in mesh.point_data)
        # the _COMP tables would have given STRESS:PLY-shaped members
        assert "STRESS:L1" not in mesh.cell_data

    def test_partial_mid_side_nodes_read_as_linear(self, engine):
        # the fixtures' CTRIA6 lists [6, 4, 60, 0, 65, 64]
        mesh = engine(path_of("static_elements"))
        eids = np.concatenate(mesh.cell_data["nastran:eid"])
        types = types_of(mesh)
        assert types[eids == 61].tolist() == ["triangle"]
        assert "triangle6" not in types

    def test_points_only_keeps_the_step_but_no_results(self, engine):
        mesh = engine(path_of("modes_elements"), points_only=True, time_step=1)
        assert mesh.field_data["nastran:mode"].tolist() == [2]
        assert not mesh.point_data
        assert sorted(mesh.cell_data) == ["nastran:eid", "nastran:pid"]

    def test_arrays_narrows_to_the_names_asked_for(self, engine):
        mesh = engine(path_of("static_elements"), arrays=["DISPLACEMENT", "STRESS:X"])
        assert list(mesh.point_data) == ["DISPLACEMENT"]
        assert sorted(mesh.cell_data) == ["STRESS:X", "nastran:eid", "nastran:pid"]

    def test_out_of_range_step_is_an_error_not_a_clamp(self, engine):
        with pytest.raises(ReadError, match="out of range"):
            engine(path_of("modes_elements"), time_step=3)
        assert engine(path_of("modes_elements"), time_step=-3).field_data[
            "nastran:mode"
        ].tolist() == [1]


def types_of(mesh):
    return np.concatenate([[c.type] * len(c.data) for c in mesh.cells])


class TestAgainstPyNastran:
    """pyNastran's own reading of the same files, frozen in ``pynastran_reference.npz``."""

    def test_every_nodal_vector_matches(self, engine):
        ref = np.load(FIXTURES / "pynastran_reference.npz")
        keys = sorted({k.rsplit("/", 1)[0] for k in ref.files})
        assert len(keys) == 38
        cache = {}
        for key in keys:
            name, table, domain = key.split("/")
            path = path_of(name)
            if name not in cache:
                with h5py.File(path) as f:
                    ids = f["/NASTRAN/INPUT/NODE/GRID"]["ID"][()]
                steps = {}
                for s in range(_sequence.num_steps(path)):
                    mesh = engine(path, time_step=s)
                    steps[int(mesh.field_data["nastran:domain"][0])] = mesh
                cache[name] = ({int(g): k for k, g in enumerate(ids)}, steps)
            pos, steps = cache[name]
            mesh = steps[int(domain)]
            rows = [
                (r, pos[int(i)])
                for r, i in enumerate(ref[key + "/ids"])
                if int(i) in pos
            ]
            r = np.array([a for a, _ in rows])
            p = np.array([b for _, b in rows])
            got = np.hstack(
                [mesh.point_data[table][p], mesh.point_data[table + "_ROT"][p]]
            )
            np.testing.assert_array_equal(got, ref[key + "/values"][r], err_msg=key)


class TestEnginesAgree:
    @pytest.mark.parametrize("name", NAMES)
    def test_every_fixture_every_step(self, name):
        path = path_of(name)
        for step in range(_sequence.num_steps(path)):
            a = _core.nastran_h5_read(path, False, None, step)
            b = py_nh5.read(path, time_step=step)
            np.testing.assert_array_equal(a.points, b.points)
            assert [(c.type, c.data.tolist()) for c in a.cells] == [
                (c.type, c.data.tolist()) for c in b.cells
            ]
            assert sorted(a.point_data) == sorted(b.point_data)
            for key in a.point_data:
                np.testing.assert_array_equal(a.point_data[key], b.point_data[key])
            assert sorted(a.cell_data) == sorted(b.cell_data)
            for key in a.cell_data:
                for x, y in zip(a.cell_data[key], b.cell_data[key]):
                    np.testing.assert_array_equal(x, y)
            assert {k: v.tolist() for k, v in a.field_data.items()} == {
                k: v.tolist() for k, v in b.field_data.items()
            }
            assert [(r.name, r.dim, r.tag, r.entries.tolist()) for r in a.regions] == [
                (r.name, r.dim, r.tag, r.entries.tolist()) for r in b.regions
            ]

    def test_the_core_reads_it_without_falling_back(self):
        set_strict_core(True)
        try:
            mesh = meshioplusplus.nastran_h5.read(
                path_of("modes_elements"), time_step=1
            )
        finally:
            set_strict_core(None)
        assert mesh.field_data["nastran:mode"].tolist() == [2]


# ---------------------------------------------------------------------------
# Built files
# ---------------------------------------------------------------------------


class TestBuiltFiles:
    def test_quadratic_mid_nodes_sit_on_their_edges(self, engine, tmp_path):
        mesh = engine(quadratic_model(tmp_path))
        assert [c.type for c in mesh.cells] == [
            "hexahedron20",
            "wedge15",
            "pyramid13",
            "quad8",
            "tetra10",
            "triangle6",
        ]
        for block in mesh.cells:
            conn = block.data[0]
            ncorner = len(conn) - len(MIDS[block.type])
            for k, (a, b) in enumerate(MIDS[block.type]):
                np.testing.assert_allclose(
                    mesh.points[conn[ncorner + k]],
                    (mesh.points[conn[a]] + mesh.points[conn[b]]) / 2,
                    err_msg=f"{block.type} mid-node {k}",
                )
            if block.type[0] in "hwpt" and block.type != "triangle6":
                assert (signed_volumes(mesh, block) > 0).all()

    def test_results_per_domain_with_nan_where_absent(self, engine, tmp_path):
        path = hex_with_results(tmp_path)
        first = engine(path)
        assert first.field_data["meshio:time"].tolist() == [0.5]
        np.testing.assert_array_equal(
            first.point_data["DISPLACEMENT"][:, 0], [1, 2, 3, 4, 5, 6, 7, 8]
        )
        second = engine(path, time_step=1)
        d = second.point_data["DISPLACEMENT"]
        assert np.isnan(d[2]).all() and d[3, 0] == -4
        assert second.point_data["DISPLACEMENT_ROT"][0, 2] == 1
        assert py_nh5.time_values(path) == [0.5, 1.5]
        assert meshioplusplus.read_metadata(path)["time_values"] == [0.5, 1.5]

    def test_read_sequence_walks_the_domains(self, tmp_path):
        steps = list(meshioplusplus.read_sequence(hex_with_results(tmp_path)))
        assert [t for t, _ in steps] == [0.5, 1.5]

    def test_a_file_without_results_has_one_step(self, engine, tmp_path):
        path = write_msc(
            tmp_path / "mesh.h5", UNIT_HEX, [("CHEXA", 20, [(1, 1, list(UNIT_HEX))])]
        )
        mesh = engine(path, time_step=-1)
        assert [c.type for c in mesh.cells] == ["hexahedron"]
        assert not mesh.field_data
        with pytest.raises(ReadError):
            engine(path, time_step=1)
        assert _sequence.num_steps(path) == 1

    def test_an_undefined_system_keeps_coordinates_and_flags_them(
        self, engine, tmp_path
    ):
        path = write_msc(tmp_path / "cp.h5", UNIT_HEX, cp={2: 5})
        mesh = engine(path)
        np.testing.assert_array_equal(mesh.points[1], [1, 0, 0])
        assert mesh.point_data["nastran:cp"].tolist() == [0, 5, 0, 0, 0, 0, 0, 0]
        assert "nastran:cd" not in mesh.point_data

    def test_elements_on_scalar_points_are_skipped(self, engine, tmp_path):
        path = write_msc(
            tmp_path / "spoint.h5",
            UNIT_HEX,
            [("CROD", 2, [(1, 1, [1, 2]), (2, 1, [1, 500])])],
            spoints=[500],
        )
        mesh = engine(path)
        assert mesh.cell_data["nastran:eid"][0].tolist() == [1]

    def test_an_undefined_grid_is_an_error(self, engine, tmp_path):
        path = write_msc(
            tmp_path / "bad.h5", UNIT_HEX, [("CROD", 2, [(1, 1, [1, 99])])]
        )
        with pytest.raises(ReadError, match="undefined GRID 99"):
            engine(path)

    def test_other_vendors_are_refused(self, engine, tmp_path):
        path = write_msc(tmp_path / "nx.h5", UNIT_HEX, version=b"NX Nastran 2019")
        with pytest.raises(ReadError, match="not MSC Nastran"):
            engine(path)

    def test_a_file_without_version_is_read(self, engine, tmp_path):
        assert (
            len(engine(write_msc(tmp_path / "nov.h5", UNIT_HEX, version=None)).points)
            == 8
        )

    def test_an_hdf5_file_that_is_not_nastran_is_refused(self, engine, tmp_path):
        path = tmp_path / "other.h5"
        with h5py.File(path, "w") as f:
            f["data"] = np.arange(3)
        with pytest.raises(ReadError, match="not an MSC Nastran HDF5"):
            engine(str(path))

    def test_a_text_file_is_refused(self, engine, tmp_path):
        path = tmp_path / "text.h5"
        path.write_text("hello\n")
        with pytest.raises(ReadError, match="not an HDF5 file"):
            engine(str(path))

    def test_a_table_without_index_is_skipped(self, engine, tmp_path):
        path = hex_with_results(tmp_path, index=False)
        mesh = engine(path)
        assert not mesh.point_data and not mesh.field_data

    def test_results_without_domains_are_an_error(self, engine, tmp_path):
        path = hex_with_results(tmp_path)
        with h5py.File(path, "a") as f:
            del f["/NASTRAN/RESULT/DOMAINS"]
        with pytest.raises(ReadError, match="DOMAINS"):
            engine(path)

    def test_an_index_naming_an_unknown_domain_is_an_error(self, engine, tmp_path):
        path = hex_with_results(tmp_path)
        with h5py.File(path, "a") as f:
            idx = f["/INDEX/NASTRAN/RESULT/NODAL/DISPLACEMENT"]
            row = idx[1]
            row["DOMAIN_ID"] = 9
            idx[1] = row
        with pytest.raises(ReadError, match="domain 9"):
            engine(path)


class TestRegistration:
    def test_read_by_extension(self):
        mesh = meshioplusplus.read(path_of("static_elements"))
        assert "DISPLACEMENT" in mesh.point_data

    def test_gid_keeps_post_h5(self):
        formats = meshioplusplus.formats()
        assert formats["extensions"][".h5"] == ["nastran_h5"]
        assert "gid" in formats["extensions"][".post.h5"]

    def test_it_is_read_only(self, tmp_path):
        formats = meshioplusplus.formats()
        assert "nastran_h5" in formats["readable"]
        assert "nastran_h5" not in formats["writable"]
        with pytest.raises(WriteError):
            meshioplusplus.write(
                tmp_path / "out.h5", meshioplusplus.read(path_of("static_elements"))
            )

    def test_core_has_the_function(self):
        assert hasattr(_core, "nastran_h5_read")


# --- coordinate systems ---------------------------------------------------------------


def _cord_probe(tmp_path):
    """``static_elements.h5`` with the GRIDs and systems of ``cord_reference.npz``
    (see tools/gen_nastran_cord_reference.py): a tilted CORD2C, a CORD2S defined
    in it and a CORD1R through three GRIDs, one of them in the CORD2C."""
    ref = np.load(FIXTURES / "cord_reference.npz")
    path = tmp_path / "cord_probe.h5"
    path.write_bytes(pathlib.Path(path_of("static_elements")).read_bytes())
    with h5py.File(path, "r+") as f:
        grid = f["NASTRAN/INPUT/NODE/GRID"][()]
        grid["X"] = ref["x_local"]
        grid["CP"] = ref["cp"]
        grid["CD"] = ref["cd"]
        del f["NASTRAN/INPUT/NODE/GRID"]
        f["NASTRAN/INPUT/NODE/GRID"] = grid
        cs = f["NASTRAN/INPUT/COORDINATE_SYSTEM"]
        dtype = cs["CORD2R"].dtype
        for name in list(cs):
            del cs[name]
        for row, ctype in zip(ref["cord2"], ref["cord2_type"]):
            table = np.zeros(1, dtype)
            table[0] = (int(row[0]), int(row[1]), *row[2:].tolist(), 1)
            cs["CORD2" + "RCS"[int(ctype) - 1]] = table
        c1 = np.zeros(
            1,
            [
                ("CID", "<i8"),
                ("G1", "<i8"),
                ("G2", "<i8"),
                ("G3", "<i8"),
                ("DOMAIN_ID", "<i8"),
            ],
        )
        c1[0] = (*ref["cord1r"][0].tolist(), 1)
        cs["CORD1R"] = c1
    return path, ref


def test_coordinate_systems_match_pynastran(engine, tmp_path):
    path, ref = _cord_probe(tmp_path)
    mesh = engine(path)
    np.testing.assert_allclose(mesh.points, ref["xyz_basic"], rtol=0, atol=1e-13)
    np.testing.assert_array_equal(mesh.point_data["nastran:cp"], ref["cp"])
    np.testing.assert_array_equal(mesh.point_data["nastran:cd"], ref["cd"])
    disp = np.hstack(
        [mesh.point_data["DISPLACEMENT"], mesh.point_data["DISPLACEMENT_ROT"]]
    )
    np.testing.assert_allclose(disp, ref["displacement_basic"], rtol=0, atol=1e-15)


def test_coordinate_systems_engines_agree(tmp_path):
    path, _ = _cord_probe(tmp_path)
    a = _core.nastran_h5_read(str(path))
    b = py_nh5.read(path)
    np.testing.assert_array_equal(a.points, b.points)
    for name in ("DISPLACEMENT", "DISPLACEMENT_ROT", "SPC_FORCE", "APPLIED_LOAD"):
        np.testing.assert_array_equal(a.point_data[name], b.point_data[name])
