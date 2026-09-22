"""CalculiX ``.frd`` results (read-only): both engines against real ``ccx`` output.

The single-element files and the multi-increment ones under ``meshes/frd`` were
written by ccx 2.23 from the decks in ``meshes/frd/decks`` (``tools/gen_frd_fixtures.py``
regenerates them); ``cgx_short.frd`` is hand written because ccx never writes the short
layout nor the cgx shell types 7-10.
"""

import io
import pathlib
import re

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core, _sequence
from meshioplusplus._exceptions import ReadError, WriteError
from meshioplusplus._fallback import set_strict_core
from meshioplusplus.frd import _frd as py_frd

FIXTURES = pathlib.Path(__file__).parent / "meshes" / "frd"

# Edge (a, b) of each mid-node, in meshio++ order, for a straight-sided element.
MIDS = {
    "hexahedron20": [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ],
    "wedge15": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
    "tetra10": [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
    "triangle6": [(0, 1), (1, 2), (2, 0)],
    "quad8": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "line3": [(0, 1)],
}
CORNERS = {
    "hexahedron20": 8,
    "wedge15": 6,
    "tetra10": 4,
    "triangle6": 3,
    "quad8": 4,
    "line3": 2,
}

# ccx element -> the cell type of the (possibly expanded) block it writes.
SINGLE = {
    "c3d8": "hexahedron",
    "c3d6": "wedge",
    "c3d4": "tetra",
    "c3d20": "hexahedron20",
    "c3d15": "wedge15",
    "c3d10": "tetra10",
    "t3d2": "line",
    "t3d3": "line3",
    # ccx expands shells and beams into solids
    "s3": "wedge",
    "s4": "hexahedron",
    "s6": "wedge15",
    "s8": "hexahedron20",
    "b31": "hexahedron",
    "b32": "hexahedron20",
}
NOT_EXPANDED = [k for k in SINGLE if k[0] in "ct"]
# ccx expands a quadratic beam into a he20 whose "mid-edge" nodes are the far end of the
# beam (its own convention, unusable as geometry), so only the other 13 are checked.
GEOMETRIC = [k for k in SINGLE if k != "b32"]


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read signature."""
    if request.param == "core":
        # the C++ reader itself, not the wrapper that falls back to Python
        def read(path, points_only=False, arrays=None, time_step=0, derived=False):
            return _core.frd_read(str(path), points_only, arrays, time_step, derived)

        return read
    return py_frd.read


def path_of(name):
    return str(FIXTURES / f"{name}.frd")


def canon(mesh):
    return (
        mesh.points.tolist(),
        [(c.type, np.asarray(c.data).tolist()) for c in mesh.cells],
        {k: np.asarray(v).tolist() for k, v in sorted(mesh.field_data.items())},
    )


def core_read(path, time_step=0, derived=False):
    """The C++ reader itself, with no fallback to hide a broken native path."""
    return _core.frd_read(str(path), False, None, time_step, derived)


def print_table(name):
    """(node, ux, uy, uz) rows of the ``.dat`` print ccx wrote for the same run."""
    text = (FIXTURES / f"{name}.dat").read_text()
    rows = re.findall(r"^\s+(\d+)((?:\s+[-+0-9.E]+){3})\s*$", text, re.M)
    return [(int(n), [float(v) for v in vals.split()]) for n, vals in rows]


class TestSingleElements:
    @pytest.mark.parametrize("name", sorted(SINGLE))
    def test_reads_the_expected_block(self, engine, name):
        mesh = engine(path_of(name))
        assert [c.type for c in mesh.cells] == [SINGLE[name]]
        assert len(mesh.cells[0].data) == 1
        for key in ("DISP", "STRESS", "TOSTRAIN", "ERROR"):
            assert key in mesh.point_data
        n = len(mesh.points)
        assert mesh.point_data["DISP"].shape == (n, 3)
        assert mesh.point_data["STRESS"].shape == (n, 6)
        assert mesh.point_data["ERROR"].shape == (n,)
        assert np.isfinite(mesh.point_data["STRESS"]).all()

    @pytest.mark.parametrize("name", sorted(GEOMETRIC))
    def test_mid_nodes_sit_on_their_edges(self, engine, name):
        """The he20/pe15/be3 permutations: a wrong one puts a mid-node off its edge."""
        mesh = engine(path_of(name))
        ctype = mesh.cells[0].type
        if ctype not in MIDS:
            pytest.skip("no mid-nodes")
        conn = mesh.cells[0].data[0]
        corners = CORNERS[ctype]
        for k, (a, b) in enumerate(MIDS[ctype]):
            mid = mesh.points[conn[corners + k]]
            expected = (mesh.points[conn[a]] + mesh.points[conn[b]]) / 2
            np.testing.assert_allclose(
                mid, expected, atol=1e-6, err_msg=f"{ctype} mid {k}"
            )

    @pytest.mark.parametrize("name", sorted(GEOMETRIC))
    def test_cells_are_not_inverted(self, engine, name):
        mesh = meshioplusplus.attach_quality(engine(path_of(name)))
        if "quality:volume" not in mesh.cell_data:
            pytest.skip("no volume")
        if mesh.cells[0].type.startswith("line"):
            pytest.skip("a line has no volume")
        assert (np.asarray(mesh.cell_data["quality:volume"][0]) > 0).all()

    @pytest.mark.parametrize("name", NOT_EXPANDED)
    def test_displacements_match_the_ccx_print(self, engine, name):
        """DISP against the 7-digit ``*NODE PRINT`` ccx wrote for the same nodes."""
        mesh = engine(path_of(name))
        table = print_table(name)
        assert table
        for node, values in table:
            np.testing.assert_allclose(
                mesh.point_data["DISP"][node - 1], values, rtol=1e-5, atol=1e-7
            )

    def test_glued_negative_values_are_split_by_column(self, engine):
        # `7-1.18144E-06` in the c3d20 file is node 7's x displacement, not a parse error
        text = (FIXTURES / "c3d20.frd").read_text()
        assert re.search(r"[0-9]-[0-9]\.\d{5}E", text)
        mesh = engine(path_of("c3d20"))
        assert mesh.point_data["DISP"][6][0] == pytest.approx(-1.18144e-06)

    def test_he20_connectivity_wraps_over_two_lines(self, engine):
        text = (FIXTURES / "c3d20.frd").read_text()
        assert (
            len(re.findall(r"^ -2 ", text.split("    3C")[1].split("1PSTEP")[0], re.M))
            == 2
        )
        assert sorted(engine(path_of("c3d20")).cells[0].data[0].tolist()) == list(
            range(20)
        )

    def test_tensor_components_are_named_in_the_documented_order(self):
        text = (FIXTURES / "c3d8.frd").read_text()
        block = text.split(" -4  STRESS")[1].split(" -1")[0]
        assert re.findall(r"^ -5  (\w+)", block, re.M) == [
            "SXX",
            "SYY",
            "SZZ",
            "SXY",
            "SYZ",
            "SZX",
        ]

    def test_uniaxial_truss_stress_is_all_xx(self, engine):
        stress = engine(path_of("t3d2")).point_data["STRESS"]
        assert (np.abs(stress[:, 0]) > 1e3).all()
        assert np.abs(stress[:, 1:]).max() < 1e-6

    def test_group_and_material_are_cell_data(self, engine):
        mesh = engine(path_of("c3d8"))
        assert mesh.cell_data["frd:group"][0].tolist() == [0]
        assert mesh.cell_data["frd:material"][0].tolist() == [1]


class TestExpansion:
    def test_shells_and_beams_are_solids_with_more_nodes(self, engine):
        # a 4-node shell is 4 nodes in the deck and an 8-node brick in the results
        assert len(engine(path_of("s4")).points) == 8
        assert len(engine(path_of("b31")).points) == 8
        assert len(engine(path_of("s8")).points) == 20
        assert len(engine(path_of("s3")).points) == 6


class TestSteps:
    def test_two_static_steps(self, engine):
        first = engine(path_of("mixed"))
        last = engine(path_of("mixed"), time_step=-1)
        assert first.field_data["meshio:time"].tolist() == [1.0]
        assert last.field_data["meshio:time"].tolist() == [2.0]
        assert first.field_data["frd:step"].tolist() == [1]
        assert last.field_data["frd:step"].tolist() == [2]
        assert first.field_data["frd:analysis"].tolist() == [0]
        # step 2 asked for no TOSTRAIN
        assert "TOSTRAIN" in first.point_data and "TOSTRAIN" not in last.point_data
        assert not np.allclose(first.point_data["DISP"], last.point_data["DISP"])

    def test_mixed_blocks_follow_the_element_order(self, engine):
        mesh = engine(path_of("mixed"))
        assert [(c.type, len(c.data)) for c in mesh.cells] == [
            ("hexahedron", 1),
            ("wedge", 1),
        ]
        assert [a.tolist() for a in mesh.cell_data["frd:group"]] == [[0], [0]]

    def test_frequency_modes_that_share_a_value_are_still_distinct_steps(self, engine):
        modes = [engine(path_of("freq"), time_step=i) for i in range(3)]
        values = [m.field_data["meshio:time"][0] for m in modes]
        assert values[0] == values[1] != values[2]
        assert [m.field_data["frd:step"][0] for m in modes] == [1, 2, 3]
        assert all(m.field_data["frd:analysis"][0] == 2 for m in modes)
        assert not np.allclose(modes[0].point_data["DISP"], modes[1].point_data["DISP"])

    def test_heat_transfer(self, engine):
        steps = [engine(path_of("heat"), time_step=i) for i in range(4)]
        assert [s.field_data["meshio:time"][0] for s in steps] == [0.5, 1.0, 1.75, 2.0]
        assert all(s.field_data["frd:analysis"][0] == 1 for s in steps)
        temp = steps[-1].point_data["NDTEMP"]
        assert temp.shape == (12,) and temp[0] == pytest.approx(100.0)
        assert steps[0].point_data["FLUX"].shape == (12, 3)

    def test_out_of_range_step_is_an_error_not_a_clamp(self, engine):
        for bad in (2, -3, 99):
            with pytest.raises(ReadError, match="out of range"):
                engine(path_of("mixed"), time_step=bad)

    def test_a_file_without_results_has_no_steps(self, engine, tmp_path):
        text = (FIXTURES / "c3d8.frd").read_text().split("    1PSTEP")[0] + "  9999\n"
        path = tmp_path / "geometry.frd"
        path.write_text(text)
        mesh = engine(str(path))
        assert not mesh.point_data and "meshio:time" not in mesh.field_data
        assert engine(str(path), time_step=-1).cells[0].type == "hexahedron"
        with pytest.raises(ReadError, match="no time steps"):
            engine(str(path), time_step=1)

    def test_metadata_lists_every_increment(self):
        meta = meshioplusplus.read_metadata(path_of("freq"))
        assert meta["format"] == "frd"
        assert len(meta["time_values"]) == 3
        assert meshioplusplus.read_metadata(path_of("mixed"))["time_values"] == [
            1.0,
            2.0,
        ]
        assert _sequence.num_steps(path_of("heat")) == 4

    def test_python_reader_attaches_the_series(self):
        assert py_frd.read(path_of("heat")).time_values == [0.5, 1.0, 1.75, 2.0]

    def test_read_sequence_walks_the_increments(self):
        steps = list(meshioplusplus.read_sequence(path_of("heat")))
        assert [t for t, _ in steps] == [0.5, 1.0, 1.75, 2.0]
        assert steps[3][1].point_data["NDTEMP"][0] == pytest.approx(100.0)


class TestShortFormat:
    """`cgx_short.frd`: I5 ids, cgx shell types 7-10, an 8-component result, NaN fill."""

    def test_shell_types_and_their_group_and_material(self, engine):
        mesh = engine(path_of("cgx_short"))
        assert [(c.type, c.data.shape) for c in mesh.cells] == [
            ("quad", (1, 4)),
            ("quad8", (1, 8)),
            ("triangle", (1, 3)),
            ("triangle6", (1, 6)),
        ]
        assert mesh.cell_data["frd:group"][0].tolist() == [2]
        assert mesh.cell_data["frd:material"][0].tolist() == [3]

    def test_shell_mid_nodes_are_where_meshio_expects_them(self, engine):
        mesh = engine(path_of("cgx_short"))
        for block in mesh.cells:
            if block.type not in MIDS:
                continue
            conn = block.data[0]
            for k, (a, b) in enumerate(MIDS[block.type]):
                np.testing.assert_allclose(
                    mesh.points[conn[CORNERS[block.type] + k]],
                    (mesh.points[conn[a]] + mesh.points[conn[b]]) / 2,
                    atol=1e-9,
                )

    def test_continuation_lines_join_into_one_row(self, engine):
        wide = engine(path_of("cgx_short")).point_data["WIDE"]
        assert wide.shape == (21, 8)
        np.testing.assert_allclose(wide[0], np.arange(8) + 10.0 - 3.5 + 0.0 * 1)
        np.testing.assert_allclose(wide[1], np.arange(8) + 20.0 - 3.5)

    def test_calculated_components_are_not_data(self, engine):
        # DISP has four -5 lines; the fourth (ALL) is iexist 1: computed, not stored
        assert engine(path_of("cgx_short")).point_data["DISP"].shape == (21, 3)

    def test_nodes_a_block_omits_are_nan(self, engine):
        temp = engine(path_of("cgx_short"), time_step=1).point_data["NDTEMP"]
        assert temp.shape == (21,)
        np.testing.assert_allclose(temp[:6], 21.0 + np.arange(6))
        assert np.isnan(temp[6:]).all()


class TestSelectiveRead:
    def test_points_only_skips_the_data_but_keeps_the_step(self, engine):
        mesh = engine(path_of("mixed"), points_only=True, time_step=1)
        assert not mesh.point_data
        assert mesh.field_data["meshio:time"].tolist() == [2.0]

    def test_arrays_narrows_to_the_names_asked_for(self, engine):
        mesh = engine(path_of("c3d8"), arrays=["DISP"])
        assert sorted(mesh.point_data) == ["DISP"]
        assert not engine(path_of("c3d8"), arrays=[]).point_data

    def test_generic_read_reaches_the_reader(self):
        mesh = meshioplusplus.read(path_of("mixed"), time_step=-1, arrays=["DISP"])
        assert sorted(mesh.point_data) == ["DISP"]
        last = meshioplusplus.read(path_of("mixed"), time_step=-1)
        assert last.field_data["meshio:time"].tolist() == [2.0]
        np.testing.assert_array_equal(last.point_data["DISP"], mesh.point_data["DISP"])


class TestDerived:
    def independent(self, s):
        s = np.asarray(s)
        mises = np.sqrt(
            0.5
            * (
                (s[:, 0] - s[:, 1]) ** 2
                + (s[:, 1] - s[:, 2]) ** 2
                + (s[:, 2] - s[:, 0]) ** 2
                + 6 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2)
            )
        )
        full = np.empty((len(s), 3, 3))
        for i, (xx, yy, zz, xy, yz, zx) in enumerate(s):
            full[i] = [[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]]
        return mises, np.linalg.eigvalsh(full)

    @pytest.mark.parametrize("name", ["c3d20", "c3d15", "s8", "t3d3"])
    def test_stress_and_strain_invariants(self, engine, name):
        mesh = engine(path_of(name), derived=True)
        for field in ("STRESS", "TOSTRAIN"):
            mises, principal = self.independent(mesh.point_data[field])
            scale = max(np.abs(mesh.point_data[field]).max(), 1e-30)
            np.testing.assert_allclose(
                mesh.point_data[f"{field}_mises"], mises, rtol=1e-9, atol=1e-9 * scale
            )
            np.testing.assert_allclose(
                mesh.point_data[f"{field}_principal"],
                principal,
                rtol=1e-9,
                atol=1e-9 * scale,
            )
            assert (np.diff(mesh.point_data[f"{field}_principal"], axis=1) >= 0).all()

    def test_off_by_default_and_only_for_tensors(self, engine):
        assert not [
            k
            for k in engine(path_of("c3d8")).point_data
            if k.endswith(("_mises", "_principal"))
        ]
        derived = engine(path_of("c3d8"), derived=True).point_data
        assert "DISP_mises" not in derived and "ERROR_mises" not in derived

    def test_arrays_can_ask_for_a_derived_name_alone(self, engine):
        mesh = engine(path_of("c3d8"), derived=True, arrays=["STRESS_mises"])
        assert sorted(mesh.point_data) == ["STRESS_mises"]

    def test_nan_nodes_stay_nan(self, engine, tmp_path):
        text = (FIXTURES / "c3d8.frd").read_text()
        # drop the last STRESS row, so that node's tensor is NaN
        head, rest = text.split(" -4  STRESS")
        lines = rest.split("\n")
        block_end = lines.index(" -3")
        rows = [i for i, ln in enumerate(lines[:block_end]) if ln.startswith(" -1")]
        del lines[rows[-1]]
        path = tmp_path / "gap.frd"
        path.write_text(head + " -4  STRESS" + "\n".join(lines))
        mesh = engine(str(path), derived=True)
        assert np.isnan(mesh.point_data["STRESS_mises"]).sum() == 1
        assert np.isnan(mesh.point_data["STRESS_principal"]).all(axis=1).sum() == 1

    def test_engines_agree(self):
        for name in ("c3d20", "s6", "mixed"):
            a = core_read(path_of(name), derived=True)
            b = py_frd.read(path_of(name), derived=True)
            assert sorted(a.point_data) == sorted(b.point_data)
            for key in a.point_data:
                np.testing.assert_allclose(
                    a.point_data[key], b.point_data[key], rtol=1e-9, atol=1e-6
                )


class TestEnginesAgree:
    @pytest.mark.parametrize("name", sorted(p.stem for p in FIXTURES.glob("*.frd")))
    def test_every_fixture_every_step(self, name):
        steps = meshioplusplus.read_metadata(path_of(name))["time_values"]
        for step in range(max(len(steps), 1)):
            a = core_read(path_of(name), time_step=step)
            b = py_frd.read(path_of(name), time_step=step)
            assert canon(a) == canon(b)
            assert sorted(a.point_data) == sorted(b.point_data)
            for key in a.point_data:
                assert a.point_data[key].shape == b.point_data[key].shape
                np.testing.assert_array_equal(a.point_data[key], b.point_data[key])
            assert sorted(a.cell_data) == sorted(b.cell_data)

    def test_the_core_reads_it_without_falling_back(self):
        set_strict_core(True)
        try:
            mesh = meshioplusplus.frd.read(path_of("mixed"), time_step=1)
        finally:
            set_strict_core(None)
        assert mesh.field_data["frd:step"].tolist() == [2]

    def test_core_has_the_function(self):
        assert hasattr(_core, "frd_read")


class TestRegistration:
    def test_read_by_extension_and_by_content(self, tmp_path):
        assert meshioplusplus.read(path_of("c3d8")).cells[0].type == "hexahedron"
        renamed = tmp_path / "results.dat.bak"
        renamed.write_bytes((FIXTURES / "c3d8.frd").read_bytes())
        assert meshioplusplus.read(renamed).cells[0].type == "hexahedron"
        assert meshioplusplus._sniff._sniff_format_py(renamed) == "frd"

    def test_it_is_read_only(self):
        formats = meshioplusplus.formats()
        assert "frd" in formats["readable"] and "frd" not in formats["writable"]
        assert formats["extensions"][".frd"] == ["frd"]
        with pytest.raises(WriteError):
            meshioplusplus.write("out.frd", meshioplusplus.read(path_of("c3d8")))

    def test_buffer_read_uses_the_python_reader(self):
        text = (FIXTURES / "mixed.frd").read_text()
        mesh = meshioplusplus.frd.read(io.StringIO(text))
        assert canon(mesh) == canon(py_frd.read(path_of("mixed")))

    def test_sniff_needs_the_user_header_or_a_node_block(self, tmp_path):
        for head, expected in (
            ("    1C\n    1UTEXT\n", "frd"),
            ("    1C\r\n    2C   \n", "frd"),
            ("    1C\n", ""),
            ("    1C\nsomething else\n", ""),
            ("1C\n    1U\n", ""),
        ):
            path = tmp_path / "probe.bin"
            path.write_text(head)
            assert meshioplusplus._sniff._sniff_format_py(path) == expected
            assert _core.sniff_format(str(path)) == expected


class TestDatFile:
    """``meshioplusplus.frd.read_dat``: the ``*NODE PRINT``/``*EL PRINT`` tabular
    print, a companion file with no mesh in it. Python-only, not registered as a
    format (``.dat`` belongs to Tecplot)."""

    def test_nodal_displacements_match_the_frd_reading_of_the_same_run(self):
        tables = py_frd.read_dat(path_of("c3d8").replace(".frd", ".dat"))
        disp = next(t for t in tables if t["quantity"] == "displacements")
        assert disp["kind"] == "node"
        assert disp["components"] == ["vx", "vy", "vz"]
        assert "int_points" not in disp
        mesh = meshioplusplus.frd.read(path_of("c3d8"))
        frd_disp = np.asarray(mesh.point_data["DISP"])
        for node_id, row in zip(disp["ids"], disp["values"]):
            np.testing.assert_allclose(frd_disp[node_id - 1], row, atol=1e-6)

    def test_element_stresses_are_per_integration_point_in_a_different_component_order(
        self,
    ):
        tables = py_frd.read_dat(path_of("c3d8").replace(".frd", ".dat"))
        stress = next(t for t in tables if t["quantity"] == "stresses")
        assert stress["kind"] == "element"
        # .dat order is sxx syy szz sxy sxz syz -- .frd's is xx yy zz xy yz zx (the
        # last two swapped).
        assert stress["components"] == ["sxx", "syy", "szz", "sxy", "sxz", "syz"]
        assert stress["set"] == "E"
        # A single C3D8 element: 8 integration points, all on element 1.
        assert list(stress["ids"]) == [1] * 8
        assert sorted(stress["int_points"]) == list(range(1, 9))
        assert stress["values"].shape == (8, 6)
        assert np.isfinite(stress["values"]).all()
        # Per-Gauss-point values differ from the nodally-extrapolated .frd ones, but
        # both describe the same physical field: element-averaged magnitudes agree
        # to within the extrapolation/averaging error.
        mesh = meshioplusplus.frd.read(path_of("c3d8"))
        dat_mean_sxx = stress["values"][:, 0].mean()
        frd_mean_xx = np.asarray(mesh.point_data["STRESS"])[:, 0].mean()
        assert dat_mean_sxx == pytest.approx(frd_mean_xx, rel=0.05)

    def test_step_and_increment_and_time_are_attached(self):
        tables = py_frd.read_dat(path_of("c3d8").replace(".frd", ".dat"))
        for t in tables:
            assert t["step"] == 1
            assert t["increment"] == 1
            assert t["time"] == pytest.approx(1.0)

    def test_every_solid_single_element_fixture_reads(self):
        # Every SOLIDS single-element deck now emits an *EL PRINT stress section
        # (tools/gen_frd_fixtures.py); a broad sweep across element types with
        # different node/integration-point counts.
        for name in ("c3d4", "c3d6", "c3d8", "c3d10", "c3d15", "c3d20"):
            tables = py_frd.read_dat(path_of(name).replace(".frd", ".dat"))
            kinds = {t["quantity"]: t for t in tables}
            assert "displacements" in kinds
            assert "stresses" in kinds
            assert kinds["stresses"]["kind"] == "element"
            assert len(kinds["stresses"]["ids"]) > 0

    def test_is_reexported_from_the_package(self):
        assert meshioplusplus.frd.read_dat is py_frd.read_dat

    def test_not_a_dat_file_fails_cleanly(self, tmp_path):
        path = tmp_path / "junk.dat"
        path.write_text("hello\nworld\n")
        with pytest.raises(ReadError, match="no recognisable"):
            py_frd.read_dat(str(path))

    def test_missing_file_fails_cleanly(self):
        with pytest.raises(ReadError, match="could not read"):
            py_frd.read_dat("/no/such/file.dat")


def edit(tmp_path, transform, source="c3d8"):
    path = tmp_path / "edited.frd"
    path.write_text(transform((FIXTURES / f"{source}.frd").read_text()))
    return str(path)


class TestNumberSpellings:
    def test_crlf_files_read_like_lf_ones(self, engine, tmp_path):
        # a Windows checkout may already hold CRLF fixtures: normalise, then convert
        raw = (FIXTURES / "mixed.frd").read_bytes().replace(b"\r\n", b"\n")
        text = raw.replace(b"\n", b"\r\n")
        path = tmp_path / "crlf.frd"
        path.write_bytes(text)
        for step in (0, 1):
            a = engine(str(path), time_step=step)
            b = engine(path_of("mixed"), time_step=step)
            assert canon(a) == canon(b)
            for key in b.point_data:
                np.testing.assert_array_equal(a.point_data[key], b.point_data[key])

    def test_doubled_carriage_returns_are_not_blank_lines(self, engine, tmp_path):
        raw = (FIXTURES / "mixed.frd").read_bytes().replace(b"\r\n", b"\n")
        path = tmp_path / "cr.frd"
        path.write_bytes(raw.replace(b"\n", b"\r\r\n"))
        assert canon(engine(str(path))) == canon(engine(path_of("mixed")))

    def test_fortran_exponents_and_special_values(self, engine, tmp_path):
        def spell(t):
            head, rest = t.split(" -4  DISP")
            first = re.search(r"^ -1         1 ([^\n]*)$", rest, re.M).group(0)
            row = " -1         1" + " 1.50000D+01" + "         NaN" + "   -Infinity"
            return head + " -4  DISP" + rest.replace(first, row, 1)

        path = edit(tmp_path, spell)
        u = engine(path).point_data["DISP"][0]
        assert u[0] == 15.0 and np.isnan(u[1]) and u[2] == -np.inf


class TestBinaryLayout:
    """The binary layout (*NODE OUTPUT/*ELEMENT OUTPUT) ccx writes: raw little-endian
    records with the ASCII headers kept, cross-checked against the ASCII rendition of
    the exact same run. See doc/formats/frd.md."""

    @pytest.mark.parametrize("name", ["c3d8", "c3d20", "cantilever_static"])
    def test_matches_the_ascii_rendition_of_the_same_run(self, engine, name):
        ascii_mesh = engine(path_of(name))
        bin_mesh = engine(path_of(f"{name}_bin"))
        np.testing.assert_allclose(ascii_mesh.points, bin_mesh.points, atol=1e-9)
        assert len(ascii_mesh.cells) == len(bin_mesh.cells)
        for ca, cb in zip(ascii_mesh.cells, bin_mesh.cells):
            assert ca.type == cb.type
            np.testing.assert_array_equal(ca.data, cb.data)
        common = set(ascii_mesh.point_data) & set(bin_mesh.point_data)
        assert common  # DISP at least, on every fixture
        for key in common:
            # ASCII is E12.5 (~6 sig figs); binary result values are float32.
            np.testing.assert_allclose(
                ascii_mesh.point_data[key],
                bin_mesh.point_data[key],
                rtol=1e-3,
                atol=1e-6,
            )
        for key in ("frd:group", "frd:material"):
            for ba, bb in zip(ascii_mesh.cell_data[key], bin_mesh.cell_data[key]):
                np.testing.assert_array_equal(ba, bb)

    def test_derived_works_on_a_binary_file(self, engine):
        mesh = engine(path_of("c3d8_bin"), derived=True)
        ascii_mesh = engine(path_of("c3d8"), derived=True)
        assert "STRESS_mises" in mesh.point_data
        assert "STRESS_principal" in mesh.point_data
        assert np.asarray(mesh.point_data["STRESS_mises"]) == pytest.approx(
            np.asarray(ascii_mesh.point_data["STRESS_mises"]), rel=1e-3
        )

    def test_a_second_node_or_element_block_is_ignored_with_a_warning(self, engine):
        # cantilever_static_bin has exactly one 2C/3C block; this only pins that the
        # binary path shares the ASCII path's "second block ignored" wiring, not a
        # real duplicate-block fixture (ccx never writes one).
        mesh = engine(path_of("cantilever_static_bin"))
        assert mesh.points.shape[0] == 117

    def test_engines_agree(self):
        for name in ("c3d8", "c3d20", "cantilever_static"):
            core = core_read(path_of(f"{name}_bin"))
            python = py_frd.read(path_of(f"{name}_bin"))
            assert canon(core) == canon(python)


class TestErrors:
    def test_a_text_file_claiming_the_binary_flag_fails_cleanly(self, engine, tmp_path):
        # A hand-edited ASCII file with the 2C flag flipped to 2: the ASCII bytes get
        # misread as binary records and fail on the first structural check (an
        # invalid element type or an out-of-bounds record) rather than crashing.
        path = edit(
            tmp_path, lambda t: re.sub(r"(    2C[^\n]*?)1(\n)", r"\g<1>2\2", t, count=1)
        )
        with pytest.raises(ReadError):
            engine(path)

    def test_not_a_frd_file(self, engine, tmp_path):
        path = tmp_path / "junk.frd"
        path.write_text("hello\nworld\n")
        with pytest.raises(ReadError, match="no node block"):
            engine(str(path))

    def test_missing_file(self, engine, tmp_path):
        with pytest.raises(ReadError):
            engine(str(tmp_path / "absent.frd"))

    def test_element_on_an_undefined_node(self, engine, tmp_path):
        path = edit(
            tmp_path,
            lambda t: t.replace(
                "         8\n -3\n    1PSTEP", "        88\n -3\n    1PSTEP", 1
            ),
        )
        with pytest.raises(ReadError, match="undefined node 88"):
            engine(path)

    def test_element_with_the_wrong_node_count(self, engine, tmp_path):
        path = edit(
            tmp_path, lambda t: t.replace("    1    0    1\n", "    3    0    1\n", 1)
        )
        with pytest.raises(ReadError, match="needs 4 nodes"):
            engine(path)

    def test_result_on_an_undefined_node(self, engine, tmp_path):
        def bad(t):
            head, rest = t.split(" -4  DISP")
            return (
                head + " -4  DISP" + rest.replace(" -1         8 ", " -1        99 ", 1)
            )

        path = edit(tmp_path, bad)
        with pytest.raises(ReadError, match="undefined node 99"):
            engine(path)

    def test_truncated_result_line(self, engine, tmp_path):
        def cut(t):
            head, rest = t.split(" -4  STRESS")
            return (
                head
                + " -4  STRESS"
                + re.sub(r"( -1         3[^\n]{24})[^\n]*\n", r"\1\n", rest, count=1)
            )

        path = edit(tmp_path, cut)
        with pytest.raises(ReadError, match="missing values"):
            engine(path)

    def test_malformed_number(self, engine, tmp_path):
        path = edit(
            tmp_path,
            lambda t: t.replace(
                " 1.00000E+00 0.00000E+00 0.00000E+00\n -1         3",
                " 1.0000xE+00 0.00000E+00 0.00000E+00\n -1         3",
                1,
            ),
        )
        with pytest.raises(ReadError, match="invalid real field"):
            engine(path)

    def test_unknown_element_types_are_skipped_with_a_warning(
        self, engine, tmp_path, caplog
    ):
        path = edit(
            tmp_path, lambda t: t.replace("    1    0    1\n", "   99    0    1\n", 1)
        )
        with caplog.at_level("WARNING"):
            mesh = engine(path)
        assert mesh.cells == []
