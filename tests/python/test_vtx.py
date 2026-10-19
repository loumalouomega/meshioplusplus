"""DOLFINx VTX output (ADIOS2 ``.bp`` directories): real DOLFINx 0.11 runs
checked against ParaView 6.1's ``ADIOS2VTXReader`` (and, for the steps ParaView
cannot read, against the raw ADIOS2 arrays); see ``meshes/vtx/README.md``.

The engine is the core when it was built with ADIOS2, else the ``adios2``
package: one process never holds both, since the package brings its own
ADIOS2 libraries."""

import pathlib
import shutil
import sys

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._exceptions import ReadError
from meshioplusplus.vtx import _vtx

HAS_CORE = getattr(_core, "__has_adios2__", False)
if not HAS_CORE:
    pytest.importorskip("adios2")
    if sys.platform == "win32":
        # The Windows `adios2` wheel dies with a native access violation in
        # `blocks_info` at random (it cannot be caught, and it takes the whole
        # pytest process down); Linux and macOS cover the twin.
        pytest.skip("adios2 wheel crashes on Windows", allow_module_level=True)

VTX = pathlib.Path(__file__).parent / "meshes" / "vtx"
NAMES = ["heat", "elastic", "quads", "heat_np2", "complex"]

# VTK_LAGRANGE_* (type, nodes) -> meshio++ type, as the reader lowers them.
LOWERED = {
    (69, 3): "triangle",
    (70, 4): "quad",
    (71, 10): "tetra10",
    (5, 3): "triangle",
}


def _reference(name):
    return np.load(VTX / f"{name}.reference.npz")


def _raw(name):
    return np.load(VTX / f"{name}.raw.npz")


def _connectivity(mesh):
    return np.concatenate([c.data.ravel() for c in mesh.cells])


@pytest.mark.parametrize("name", NAMES)
def test_matches_paraview(name):
    ref = _reference(name)
    assert meshioplusplus.vtx.time_values(VTX / f"{name}.bp") == list(ref["times"])
    for k in ref["steps"]:
        mesh = meshioplusplus.read(VTX / f"{name}.bp", time_step=int(k))
        np.testing.assert_array_equal(mesh.points, ref[f"k{k}/points"])
        np.testing.assert_array_equal(_connectivity(mesh), ref[f"k{k}/connectivity"])
        sizes = np.diff(ref[f"k{k}/offsets"])
        types = [LOWERED[(t, n)] for t, n in zip(ref[f"k{k}/types"], sizes)]
        assert [c.type for c in mesh.cells] == sorted(set(types), key=types.index)
        assert mesh.field_data["meshio:time"][0] == ref["times"][k]
        for key in ref.files:
            if not key.startswith(f"k{k}/point/") and not key.startswith(f"k{k}/cell/"):
                continue
            kind, array = key.split("/")[1:]
            expected = ref[key]
            if kind == "point":
                got = mesh.point_data[array]
                # ParaView appends DOLFINx's stray one-value ghost blocks.
                expected = expected[: len(got)]
            else:
                got = np.concatenate(mesh.cell_data[array])
            np.testing.assert_array_equal(got, expected)


@pytest.mark.parametrize("name", ["heat", "heat_np2", "elastic", "complex"])
def test_every_step_matches_raw_arrays(name):
    """Every step, the ones ParaView aborts on included: the mesh carried from
    step 0 (`reuse`), the functions as ADIOS2 stored them."""
    raw = _raw(name)
    first = meshioplusplus.read(VTX / f"{name}.bp")
    steps = sorted({int(key.split("/")[0][1:]) for key in raw.files})
    for k in steps:
        mesh = meshioplusplus.read(VTX / f"{name}.bp", time_step=k)
        np.testing.assert_array_equal(mesh.points, first.points)
        np.testing.assert_array_equal(_connectivity(mesh), _connectivity(first))
        for key in raw.files:
            if not key.startswith(f"k{k}/"):
                continue
            array = key.split("/")[1]
            got = (
                mesh.point_data[array]
                if array in mesh.point_data
                else np.concatenate(mesh.cell_data[array])
            )
            np.testing.assert_array_equal(got, raw[key])


def test_ghost_arrays_skip_dolfinx_stray_blocks():
    """The first step of a `reuse` file holds malformed one-value ghost blocks;
    the arrays match the rank blocks and every step's are the same."""
    ids = [
        meshioplusplus.read(VTX / "heat_np2.bp", time_step=k).point_data[
            "vtkOriginalPointIds"
        ]
        for k in range(3)
    ]
    assert len(ids[0]) == 38
    for other in ids[1:]:
        np.testing.assert_array_equal(ids[0], other)
    assert set(ids[0]) == set(range(25))


def _by_position(mesh):
    return {tuple(np.round(p, 9)): i for i, p in enumerate(mesh.points)}


@pytest.mark.parametrize("step", [0, 2])
def test_drop_ghosts_gives_the_serial_mesh(step):
    parallel = meshioplusplus.vtx.read(
        VTX / "heat_np2.bp", time_step=step, ghosts="drop"
    )
    serial = meshioplusplus.read(VTX / "heat.bp", time_step=step)
    assert len(parallel.points) == len(serial.points)
    assert "vtkGhostType" not in parallel.point_data
    index = _by_position(parallel)
    order = np.array([index[tuple(np.round(p, 9))] for p in serial.points])
    np.testing.assert_allclose(
        parallel.point_data["u"][order], serial.point_data["u"], atol=1e-14
    )

    def cells(mesh):
        return {
            tuple(sorted(tuple(np.round(mesh.points[i], 9)) for i in c))
            for c in mesh.cells[0].data
        }

    assert cells(parallel) == cells(serial)


def test_points_only_and_array_filter():
    mesh = meshioplusplus.vtx.read(VTX / "elastic.bp", points_only=True)
    assert not mesh.point_data and not mesh.cell_data
    mesh = meshioplusplus.vtx.read(VTX / "elastic.bp", arrays=["density"])
    assert list(mesh.cell_data) == ["density"] and not mesh.point_data


def test_negative_and_out_of_range_steps():
    last = meshioplusplus.read(VTX / "heat.bp", time_step=-1)
    assert last.field_data["meshio:time"][0] == pytest.approx(0.02)
    with pytest.raises(ReadError, match="3"):
        meshioplusplus.read(VTX / "heat.bp", time_step=3)


def test_metadata_and_sequence():
    meta = meshioplusplus.read_metadata(VTX / "heat.bp")
    assert meta["time_values"] == pytest.approx([0.0, 0.01, 0.02])
    steps = list(meshioplusplus.read_sequence(VTX / "heat.bp"))
    assert [t for t, _ in steps] == pytest.approx([0.0, 0.01, 0.02])


def test_glob_keeps_bp_directories(tmp_path):
    for k in range(2):
        shutil.copytree(VTX / "quads.bp", tmp_path / f"run_{k}.bp")
    (tmp_path / "ignored").mkdir()
    entries = meshioplusplus.sequence_entries(str(tmp_path / "run_*.bp"))
    assert [pathlib.Path(e["path"]).name for e in entries] == ["run_0.bp", "run_1.bp"]


def test_sniffed_as_a_directory(tmp_path):
    target = tmp_path / "noext"
    shutil.copytree(VTX / "heat.bp", target)
    from meshioplusplus._sniff import sniff_format

    assert sniff_format(target) == "vtx"


def test_refuses_a_bp_without_schema():
    with pytest.raises(ReadError, match="adios4dolfinx"):
        meshioplusplus.read(VTX / "plain.bp")


def test_python_twin_names_its_dependency(monkeypatch):
    monkeypatch.setitem(sys.modules, "adios2", None)
    with pytest.raises(ImportError, match=r"meshioplusplus\[adios2\]"):
        _vtx.read(VTX / "heat.bp")


# ---------------------------------------------------------------------------
# The writer: VTXWriter's layout for one rank, from the core when it has ADIOS2
# and from the `adios2` package otherwise.
# ---------------------------------------------------------------------------

import subprocess  # noqa: E402,F401  (the cross-engine check below)
import textwrap  # noqa: E402

from meshioplusplus import WriteError  # noqa: E402


def _tri6_mesh():
    pts = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0.5, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0]],
        dtype=float,
    )
    mesh = meshioplusplus.Mesh(pts, [("triangle6", np.array([[0, 1, 2, 3, 4, 5]]))])
    mesh.point_data["u"] = np.arange(6.0)
    mesh.point_data["v"] = np.arange(18.0).reshape(6, 3)
    mesh.cell_data["rho"] = [np.array([2.5])]
    mesh.field_data["meshio:time"] = np.array([0.25])
    return mesh


@pytest.mark.parametrize(
    "ctype,corners",
    [
        ("line", [[0, 0, 0], [1, 0, 0]]),
        ("triangle", [[0, 0, 0], [1, 0, 0], [0, 1, 0]]),
        ("quad", [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]),
        ("tetra", [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        (
            "wedge",
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]],
        ),
    ],
)
def test_each_cell_type_survives_a_round_trip(ctype, corners, tmp_path):
    mesh = meshioplusplus.Mesh(
        np.array(corners, dtype=float), [(ctype, np.arange(len(corners))[None])]
    )
    mesh.point_data["p"] = np.arange(len(corners), dtype=float)
    path = tmp_path / "m.bp"
    meshioplusplus.write(path, mesh, file_format="vtx")
    again = meshioplusplus.read(path)
    np.testing.assert_array_equal(again.points, mesh.points)
    assert [(c.type, c.data.tolist()) for c in again.cells] == [
        (ctype, [list(range(len(corners)))])
    ]
    np.testing.assert_array_equal(again.point_data["p"], mesh.point_data["p"])


def test_a_quadratic_mesh_with_data_survives(tmp_path):
    mesh = _tri6_mesh()
    path = tmp_path / "t.bp"
    meshioplusplus.write(path, mesh, file_format="vtx")
    again = meshioplusplus.read(path)
    assert [(c.type, c.data.tolist()) for c in again.cells] == [
        ("triangle6", [[0, 1, 2, 3, 4, 5]])
    ]
    np.testing.assert_array_equal(again.point_data["u"], mesh.point_data["u"])
    np.testing.assert_array_equal(again.point_data["v"], mesh.point_data["v"])
    np.testing.assert_array_equal(again.cell_data["rho"][0], [2.5])
    assert again.field_data["meshio:time"][0] == 0.25


def test_the_layout_is_vtxwriters(tmp_path):
    """The attribute and the variables a reader (ParaView's included) looks for."""
    adios2 = pytest.importorskip("adios2")
    path = tmp_path / "layout.bp"
    meshioplusplus.write(path, _tri6_mesh(), file_format="vtx")
    with adios2.Stream(str(path), "r") as stream:
        for _ in stream.steps():
            assert "vtk.xml" in stream.available_attributes()
            xml = stream.read_attribute_string("vtk.xml")
            xml = "".join(xml) if isinstance(xml, (list, tuple)) else xml
            variables = stream.available_variables()
            for name in ("geometry", "connectivity", "types", "step", "u", "v", "rho"):
                assert name in variables, name
            assert variables["types"]["SingleValue"] == "true"
            assert variables["types"]["Min"] == "69"  # VTK_LAGRANGE_TRIANGLE
            assert variables["geometry"]["Shape"] == ""  # a rank-local array
            assert 'Name="geometry"' in xml and "<CellData>" in xml
            break


def test_a_series_is_a_step_per_mesh(tmp_path):
    steps = []
    for k in range(3):
        mesh = _tri6_mesh()
        mesh.point_data["u"] = mesh.point_data["u"] + 10 * k
        steps.append((0.5 * k, mesh))
    path = tmp_path / "series.bp"
    meshioplusplus.write_sequence(path, iter(steps), file_format="vtx")
    assert meshioplusplus.vtx.time_values(path) == [0.0, 0.5, 1.0]
    for k in range(3):
        again = meshioplusplus.read(path, time_step=k)
        np.testing.assert_array_equal(
            again.point_data["u"], steps[k][1].point_data["u"]
        )
        assert again.field_data["meshio:time"][0] == 0.5 * k


def test_refusals(tmp_path):
    empty = meshioplusplus.Mesh(np.zeros((2, 3)), [])
    with pytest.raises(WriteError, match="no cells"):
        meshioplusplus.write(tmp_path / "empty.bp", empty, file_format="vtx")
    other = _tri6_mesh()
    other.point_data["extra"] = np.zeros(6)
    with pytest.raises(WriteError):
        with meshioplusplus.vtx.SeriesWriter(tmp_path / "bad.bp") as w:
            w.write(0.0, _tri6_mesh())
            w.write(1.0, other)
    with pytest.raises(ValueError, match="mesh_policy"):
        meshioplusplus.vtx.SeriesWriter(tmp_path / "policy.bp", mesh_policy="never")


def _mixed_mesh(u0=0.0, x0=0.0):
    pts = np.array(
        [[x0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0], [2, 1, 0], [5, 5, 5]],
        dtype=float,
    )
    mesh = meshioplusplus.Mesh(
        pts,
        [
            ("triangle", np.array([[0, 1, 2]])),
            ("quad", np.array([[1, 4, 5, 2]])),
            ("triangle", np.array([[0, 2, 3]])),
        ],
    )
    u = np.arange(7.0)
    u[0] = u0
    mesh.point_data["u"] = u
    mesh.point_data["v"] = np.arange(21.0).reshape(7, 3)
    mesh.cell_data["rho"] = [np.array([1.0]), np.array([2.0]), np.array([3.0])]
    mesh.field_data["meshio:time"] = np.array([0.5])
    return mesh


def test_several_cell_types_are_one_block_each_and_read_back_welded(tmp_path):
    """One block per cell type, each with the points it uses; the shared points are
    ghost copies the reader welds back, in the mesh's order, the unused point too."""
    mesh = _mixed_mesh()
    path = tmp_path / "mixed.bp"
    meshioplusplus.write(path, mesh, file_format="vtx")
    back = meshioplusplus.read(path)
    np.testing.assert_array_equal(back.points, mesh.points)
    assert [(c.type, c.data.tolist()) for c in back.cells] == [
        (c.type, c.data.tolist()) for c in mesh.cells
    ]
    np.testing.assert_array_equal(back.point_data["u"], mesh.point_data["u"])
    np.testing.assert_array_equal(back.point_data["v"], mesh.point_data["v"])
    for a, b in zip(back.cell_data["rho"], mesh.cell_data["rho"]):
        np.testing.assert_array_equal(a, b)
    assert "vtkGhostType" not in back.point_data
    assert "vtkOriginalPointIds" not in back.point_data
    assert back.field_data["meshio:time"][0] == 0.5


@pytest.mark.skipif(HAS_CORE, reason="reads the raw file with the adios2 package")
def test_the_blocks_of_a_mixed_file_hold_the_points_they_use(tmp_path):
    adios2 = pytest.importorskip("adios2")
    path = tmp_path / "mixed.bp"
    meshioplusplus.write(path, _mixed_mesh(), file_format="vtx")
    with adios2.FileReader(str(path)) as reader:
        assert "meshio.shared_points" in reader.available_attributes()
        info = reader.all_blocks_info("geometry")[0]
        counts = [int(i["Count"].split(",")[0]) for i in info]
        # triangle (0 1 2), quad (1 4 5 2), triangle (0 2 3) and the unused point 6
        assert counts == [3 + 1, 4, 3]
        ghosts = reader.all_blocks_info("vtkGhostType")[0]
        assert len(ghosts) == 3
        types = [int(float(i["Min"])) for i in reader.all_blocks_info("types")[0]]
        assert types == [69, 70, 69]  # VTK_LAGRANGE_TRIANGLE, _QUADRILATERAL


def test_parts_are_one_block_each(tmp_path):
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0]], dtype=float)
    mesh = meshioplusplus.Mesh(
        pts, [("triangle", np.array([[0, 1, 2], [0, 2, 3], [1, 4, 2]]))]
    )
    mesh.cell_data["partition:part"] = [np.array([0, 1, 0])]
    mesh.cell_data["rho"] = [np.array([10.0, 20.0, 30.0])]
    path = tmp_path / "parts.bp"
    meshioplusplus.write(path, mesh, file_format="vtx")
    back = meshioplusplus.read(path)
    np.testing.assert_array_equal(back.points, pts)
    # the cells come back grouped by part: part 0's (0, 2), then part 1's (1)
    assert back.cells[0].data.tolist() == [[0, 1, 2], [1, 4, 2], [0, 2, 3]]
    np.testing.assert_array_equal(back.cell_data["partition:part"][0], [0, 0, 1])
    np.testing.assert_array_equal(back.cell_data["rho"][0], [10.0, 30.0, 20.0])


def _steps(count=3):
    steps = []
    for k in range(count):
        mesh = _mixed_mesh(u0=10.0 * k)
        mesh.field_data["meshio:time"] = np.array([0.5 * k])
        steps.append((0.5 * k, mesh))
    return steps


def test_the_reuse_policy_writes_the_mesh_once(tmp_path):
    path = tmp_path / "reuse.bp"
    steps = _steps()
    with meshioplusplus.vtx.SeriesWriter(path, mesh_policy="reuse") as writer:
        for time, mesh in steps:
            writer.write(time, mesh)
    assert meshioplusplus.vtx.time_values(path) == [0.0, 0.5, 1.0]
    for k, (_, mesh) in enumerate(steps):
        back = meshioplusplus.read(path, time_step=k)
        np.testing.assert_array_equal(back.points, mesh.points)
        assert [c.type for c in back.cells] == [c.type for c in mesh.cells]
        np.testing.assert_array_equal(back.point_data["u"], mesh.point_data["u"])
    # the policy goes through write_sequence, and the file is smaller than update's
    other = tmp_path / "update.bp"
    meshioplusplus.write_sequence(other, iter(_steps()), file_format="vtx")
    again = tmp_path / "again.bp"
    meshioplusplus.write_sequence(
        again, iter(_steps()), file_format="vtx", mesh_policy="reuse"
    )

    def size(directory):
        return sum(f.stat().st_size for f in directory.rglob("*") if f.is_file())

    assert size(again) <= size(other)


@pytest.mark.skipif(HAS_CORE, reason="reads the raw file with the adios2 package")
def test_a_reuse_file_has_the_mesh_in_step_zero_only(tmp_path):
    adios2 = pytest.importorskip("adios2")
    path = tmp_path / "reuse.bp"
    with meshioplusplus.vtx.SeriesWriter(path, mesh_policy="reuse") as writer:
        for time, mesh in _steps():
            writer.write(time, mesh)
    with adios2.Stream(str(path), "r") as stream:
        seen = []
        for _ in stream.steps():
            seen.append(set(stream.available_variables()))
    for name in ("geometry", "connectivity", "types", "vtkGhostType"):
        assert name in seen[0]
        assert name not in seen[1] and name not in seen[2]
    for name in ("step", "u", "rho"):
        assert all(name in step for step in seen)


def test_the_reuse_policy_refuses_a_moved_mesh(tmp_path):
    with pytest.raises(WriteError, match="reuse policy"):
        with meshioplusplus.vtx.SeriesWriter(
            tmp_path / "moved.bp", mesh_policy="reuse"
        ) as writer:
            writer.write(0.0, _mixed_mesh())
            writer.write(1.0, _mixed_mesh(x0=9.0))


@pytest.mark.skipif(not HAS_CORE, reason="compares the core with the Python twin")
def test_the_python_twin_writes_what_the_core_does(tmp_path):
    """The twin runs in a child process: it loads the `adios2` package's own
    ADIOS2 libraries, which this process (holding the core's) must not."""
    pytest.importorskip("adios2")
    meshioplusplus.write(tmp_path / "core.bp", _tri6_mesh(), file_format="vtx")
    script = textwrap.dedent(
        """
        import sys, numpy as np
        import meshioplusplus
        from meshioplusplus.vtx import _write
        pts = np.array([[0,0,0],[1,0,0],[0,1,0],[0.5,0,0],[0.5,0.5,0],[0,0.5,0]], float)
        mesh = meshioplusplus.Mesh(pts, [("triangle6", np.array([[0,1,2,3,4,5]]))])
        mesh.point_data["u"] = np.arange(6.0)
        mesh.point_data["v"] = np.arange(18.0).reshape(6, 3)
        mesh.cell_data["rho"] = [np.array([2.5])]
        mesh.field_data["meshio:time"] = np.array([0.25])
        _write.write(sys.argv[1], mesh)
        """
    )
    subprocess.run(
        [sys.executable, "-c", script, str(tmp_path / "twin.bp")], check=True
    )
    a = meshioplusplus.read(tmp_path / "core.bp")
    b = meshioplusplus.read(tmp_path / "twin.bp")
    np.testing.assert_array_equal(a.points, b.points)
    assert [(c.type, c.data.tolist()) for c in a.cells] == [
        (c.type, c.data.tolist()) for c in b.cells
    ]
    for name in a.point_data:
        np.testing.assert_array_equal(a.point_data[name], b.point_data[name])
    np.testing.assert_array_equal(a.cell_data["rho"][0], b.cell_data["rho"][0])


@pytest.mark.skipif(not HAS_CORE, reason="compares the core with the Python twin")
def test_the_twin_writes_the_blocks_and_the_reuse_policy_the_core_reads(tmp_path):
    """Several cell types, parts and ``mesh_policy="reuse"``: the twin (in a child
    process) writes them, the core reads them as the mesh went in."""
    pytest.importorskip("adios2")
    script = textwrap.dedent(
        """
        import sys, numpy as np
        import meshioplusplus
        from meshioplusplus.vtx import _write
        pts = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[2,0,0],[2,1,0],[5,5,5]], float)
        def make(u0):
            mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0,1,2]])),
                ("quad", np.array([[1,4,5,2]])), ("triangle", np.array([[0,2,3]]))])
            u = np.arange(7.0); u[0] = u0
            mesh.point_data["u"] = u
            mesh.cell_data["rho"] = [np.array([1.0]), np.array([2.0]), np.array([3.0])]
            mesh.field_data["meshio:time"] = np.array([0.5])
            return mesh
        _write.write(sys.argv[1] + "/mixed.bp", make(0.0))
        w = _write.SeriesWriter(sys.argv[1] + "/reuse.bp", mesh_policy="reuse")
        for k in range(3):
            w.write(0.5 * k, make(10.0 * k))
        w.close()
        """
    )
    subprocess.run([sys.executable, "-c", script, str(tmp_path)], check=True)
    mixed = meshioplusplus.read(tmp_path / "mixed.bp")
    assert len(mixed.points) == 7
    assert [c.type for c in mixed.cells] == ["triangle", "quad", "triangle"]
    assert mixed.cells[1].data.tolist() == [[1, 4, 5, 2]]
    np.testing.assert_array_equal(mixed.point_data["u"], np.arange(7.0))
    for k in range(3):
        back = meshioplusplus.read(tmp_path / "reuse.bp", time_step=k)
        assert len(back.points) == 7
        assert back.point_data["u"][0] == 10.0 * k
        assert [c.type for c in back.cells] == ["triangle", "quad", "triangle"]
    # and the other way: the core's files read by the twin, in a child process
    core = _mixed_mesh()
    meshioplusplus.write(tmp_path / "core.bp", core, file_format="vtx")
    check = textwrap.dedent(
        """
        import sys, numpy as np
        from meshioplusplus.vtx import _vtx
        mesh = _vtx.read(sys.argv[1] + "/core.bp")
        assert len(mesh.points) == 7, len(mesh.points)
        assert [c.type for c in mesh.cells] == ["triangle", "quad", "triangle"]
        assert mesh.cells[1].data.tolist() == [[1, 4, 5, 2]]
        """
    )
    subprocess.run([sys.executable, "-c", check, str(tmp_path)], check=True)
