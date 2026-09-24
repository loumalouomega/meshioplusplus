"""Ansys MAPDL results (``.rst``, ``.rth``): both engines against files MAPDL
wrote (and pymapdl-reader's frozen reading of them: nodal solutions, averaged
element stresses and strains, reaction forces, the element records read as they
are written, a distributed solve, the full rotor of static and modal cyclic
models), plus a synthetic file that pins what those
files do not exercise -- a node rotated about all three axes, a result set
holding only some nodes, MAPDL's undefined value, a rotated element, a layered
shell, an all-zero record and the refusals."""

import collections
import math
import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.ansys_rst import _ansys_rst as py_rst

RST = pathlib.Path(__file__).parent / "meshes" / "ansys" / "rst"
FIXTURES = sorted(RST.glob("*.rst")) + sorted(RST.glob("*.rth"))
DIST = RST / "dist_static"
CYCLIC = [RST / "cyc12.rst", RST / "cyclic_v182.rst"]

_VTK = {
    1: "vertex",
    3: "line",
    5: "triangle",
    9: "quad",
    10: "tetra",
    12: "hexahedron",
    13: "wedge",
    14: "pyramid",
    21: "line3",
    22: "triangle6",
    23: "quad8",
    24: "tetra10",
    25: "hexahedron20",
    26: "wedge15",
    27: "pyramid13",
}
# pymapdl-reader's DOF labels as (meshio++ array, component).
_DOF = {
    "UX": ("U", 0),
    "UY": ("U", 1),
    "UZ": ("U", 2),
    "ROTX": ("ROT", 0),
    "ROTY": ("ROT", 1),
    "ROTZ": ("ROT", 2),
}


class _PyCyclic:
    """The Python engine's full-rotor read behind ``read_cyclic``."""

    time_values = staticmethod(py_rst.time_values)

    @staticmethod
    def read(filename, **kwargs):
        return py_rst.read(filename, **kwargs)

    @staticmethod
    def read_cyclic(filename, **kwargs):
        return py_rst.read(filename, cyclic=True, **kwargs)


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read signature."""
    if request.param == "core":
        return meshioplusplus.ansys_rst
    return _PyCyclic


def _node_index(path):
    """Node number -> point index of the mesh both engines build."""
    return py_rst._Model(py_rst._open(str(path)), True).node_index


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    for data_a, data_b in ((a.point_data, b.point_data), (a.field_data, b.field_data)):
        assert sorted(data_a) == sorted(data_b)
        for name in data_a:
            np.testing.assert_array_equal(data_a[name], data_b[name])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_array_equal(x, y)
    assert sorted((r.kind, r.name, tuple(r.entries)) for r in a.regions) == sorted(
        (r.kind, r.name, tuple(r.entries)) for r in b.regions
    )


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_engines_agree_on_every_set(path):
    times = meshioplusplus.ansys_rst.time_values(path)
    assert times == py_rst.time_values(path)
    for step in range(len(times)):
        _same(
            meshioplusplus.ansys_rst.read(path, time_step=step, lenient=True),
            py_rst.read(path, time_step=step, lenient=True),
        )


def _cells(points, cell_type, rows):
    out = collections.Counter()
    for row in rows:
        out[(cell_type, tuple(sorted(tuple(np.round(p, 7)) for p in points[row])))] += 1
    return out


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_matches_pymapdl_reader(engine, path):
    """The cells, set times and nodal solutions pymapdl-reader reads."""
    ref = np.load(RST / "pymapdl_reference.npz")
    key = path.name
    times = ref[f"{key}/times"]
    mesh = engine.read(path, lenient=True)

    expected = collections.Counter()
    points, offsets = ref[f"{key}/points"], ref[f"{key}/offsets"]
    conn = ref[f"{key}/connectivity"]
    for c, vtk in enumerate(ref[f"{key}/celltypes"]):
        if vtk:  # pymapdl-reader keeps elements it cannot map as empty cells
            row = conn[offsets[c] : offsets[c + 1]]
            expected += _cells(points, _VTK[int(vtk)], [row])
    got = collections.Counter()
    for block in mesh.cells:
        got += _cells(mesh.points, block.type, block.data)
    assert got == expected

    index = _node_index(path)
    for step, time in enumerate(times):
        mesh = engine.read(path, time_step=step, lenient=True)
        assert mesh.field_data["meshio:time"][0] == time
        rows = [index[int(n)] for n in ref[f"{key}/{step}/nnum"]]
        values = ref[f"{key}/{step}/values"]
        for k, label in enumerate(ref[f"{key}/{step}/dofs"].tolist()):
            name, component = _DOF.get(label, (label, None))
            ours = mesh.point_data[name][rows]
            if component is not None:
                ours = ours[:, component]
            np.testing.assert_allclose(ours, values[:, k], rtol=1e-12, atol=0)


# pymapdl-reader's labels of reaction DOFs as (meshio++ array, component).
_REACTION = {
    "UX": ("RF", 0),
    "UY": ("RF", 1),
    "UZ": ("RF", 2),
    "ROTX": ("RMOM", 0),
    "ROTY": ("RMOM", 1),
    "ROTZ": ("RMOM", 2),
}


def _element_keys():
    ref = np.load(RST / "pymapdl_reference.npz")
    keys = []
    for name in ref.files:
        key, _, rest = name.rpartition("/")
        if rest in ("stress", "strain", "rf_values"):
            keys.append((key.rsplit("/", 1)[0], int(key.rsplit("/", 1)[1]), rest))
    return sorted(keys)


@pytest.mark.parametrize(
    "key, step, what",
    _element_keys(),
    ids=["/".join(map(str, k)) for k in _element_keys()],
)
def test_element_results_match_pymapdl_reader(engine, key, step, what):
    """Node-averaged stresses and elastic strains (``nodal_stress``,
    ``nodal_elastic_strain``) and reaction forces (``nodal_reaction_forces``,
    which pymapdl-reader leaves in the nodal coordinate systems)."""
    ref = np.load(RST / "pymapdl_reference.npz")
    path = RST / key
    if what == "stress" and key == "temp_v13.rst":
        pytest.skip("pymapdl-reader does not read element results before v14.5")
    if what != "rf_values" and key == "beam44.rst":
        pytest.skip("line elements carry no element nodal stresses in meshio++")
    mesh = engine.read(path, time_step=step, lenient=True)
    index = _node_index(path)
    if what == "rf_values":
        model = py_rst._Model(py_rst._open(str(path)), True)
        nodal = {}
        for label, number, value in zip(
            ref[f"{key}/{step}/rf_dofs"].tolist(),
            ref[f"{key}/{step}/rf_nnum"],
            ref[f"{key}/{step}/rf_values"],
        ):
            name, component = _REACTION.get(label, (f"RF_{label}", None))
            point = index[int(number)]
            if component is None:
                assert mesh.point_data[name][point] == value
                continue
            nodal.setdefault((name, point), np.zeros(3))[component] = value
        for (name, point), vec in nodal.items():
            vec = vec.reshape(1, 3)
            py_rst._rotate(vec, model.angles[[point]])
            np.testing.assert_allclose(
                mesh.point_data[name][point], vec[0], rtol=1e-12, atol=1e-12
            )
        return
    name = "S" if what == "stress" else "EPEL"
    rows = [index[int(n)] for n in ref[f"{key}/{step}/{what}_nnum"]]
    expected = ref[f"{key}/{step}/{what}"]
    np.testing.assert_allclose(
        mesh.point_data[name][rows], expected, rtol=1e-12, atol=0
    )


def test_v13_stress_records_are_consistent():
    """A release 13 file stores 11 items per node (``SX .. SXZ, S1 S2 S3, SINT,
    SEQV``): the principal stresses are the eigenvalues of the first six."""
    results = py_rst._Results(str(RST / "temp_v13.rst"))
    results.deck()
    assert not results.sparse_ens
    base, offsets = results.element_tables(0)
    checked = 0
    for off in offsets[:20]:
        values = results.element_record(base + int(off), py_rst._ENS)
        for row in values.reshape(-1, 11):
            xx, yy, zz, xy, yz, xz = row[:6]
            t = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
            principal = np.sort(np.linalg.eigvalsh(t))[::-1]
            np.testing.assert_allclose(
                principal, row[6:9], rtol=1e-5, atol=1e-6 * abs(row).max()
            )
            assert row[9] == pytest.approx(row[6] - row[8], rel=1e-5)
            checked += 1
    assert checked == 160


def test_static_solid186_stresses(engine):
    """The done-when of roadmap item 1.2: a static SOLID186 result's nodal
    stresses match pymapdl-reader's ``nodal_stress``; per element node they are
    cell data, NaN at the midside nodes."""
    ref = np.load(RST / "pymapdl_reference.npz")
    path = RST / "beam_static_bc.rst"
    mesh = engine.read(path)
    index = _node_index(path)
    rows = [index[int(n)] for n in ref["beam_static_bc.rst/0/stress_nnum"]]
    np.testing.assert_array_equal(
        mesh.point_data["S"][rows], ref["beam_static_bc.rst/0/stress"]
    )
    # Per element node, flattened point-major: (cells, nodes * 6), and its
    # (nodes, components) is the layout.
    assert mesh.cell_data["S"][0].shape == (40, 120)
    np.testing.assert_array_equal(mesh.field_data["ansys:layout:S"], [20, 6])
    cells = mesh.cell_data["S"][0].reshape(40, 20, 6)
    assert np.isfinite(cells[:, :8]).all() and np.isnan(cells[:, 8:]).all()
    assert mesh.cell_data["ENF"][0].shape == (40, 60)
    np.testing.assert_array_equal(mesh.field_data["ansys:layout:ENF"], [20, 3])
    assert np.isfinite(mesh.cell_data["ENF"][0]).all()
    # Nodal equilibrium: at each node the element nodal forces sum to minus the
    # reaction, or to minus the applied load (FX 20, FY 30, FZ 40 at three
    # nodes), and to zero everywhere else.
    total = np.zeros((len(mesh.points), 3))
    np.add.at(
        total, mesh.cells[0].data.ravel(), mesh.cell_data["ENF"][0].reshape(-1, 3)
    )
    residual = total + np.nan_to_num(mesh.point_data["RF"])
    scale = np.abs(np.nan_to_num(mesh.point_data["RF"])).max()
    loaded = np.abs(residual) > 1e-6 * scale
    assert sorted(residual[loaded].round(5).tolist()) == [-40.0, -30.0, -20.0]


def test_distributed_solve_reads_its_partial_files(engine, tmp_path):
    """file0.rst reads file1..3.rst with it: the same model and results as the
    combined file MAPDL wrote (compared by node number)."""
    merged = engine.read(DIST / "file0.rst")
    combined = engine.read(DIST / "file.rst")
    assert len(merged.points) == len(combined.points) == 1309
    assert sorted(merged.point_data) == sorted(combined.point_data)
    a, b = _node_index(DIST / "file0.rst"), _node_index(DIST / "file.rst")
    numbers = sorted(b)
    rows_a, rows_b = [a[n] for n in numbers], [b[n] for n in numbers]
    np.testing.assert_array_equal(merged.points[rows_a], combined.points[rows_b])
    for name in combined.point_data:
        np.testing.assert_allclose(
            merged.point_data[name][rows_a],
            combined.point_data[name][rows_b],
            rtol=1e-6,
            atol=1e-9 * np.nanmax(np.abs(combined.point_data[name])),
            err_msg=name,
        )
    got = collections.Counter(b.type for b in merged.cells for _ in b.data)
    assert got == collections.Counter(b.type for b in combined.cells for _ in b.data)

    with pytest.raises(meshioplusplus.ReadError, match="not the main one"):
        engine.read(DIST / "file1.rst")
    for k in range(3):  # file3.rst left out
        (tmp_path / f"job{k}.rst").write_bytes((DIST / f"file{k}.rst").read_bytes())
    with pytest.raises(meshioplusplus.ReadError, match="one is missing"):
        engine.read(tmp_path / "job0.rst")
    (tmp_path / "other.rst").write_bytes((DIST / "file0.rst").read_bytes())
    with pytest.raises(meshioplusplus.ReadError, match="must be named"):
        engine.read(tmp_path / "other.rst")


@pytest.mark.parametrize("path", CYCLIC, ids=[p.name for p in CYCLIC])
def test_cyclic_full_rotor_matches_pymapdl_reader(engine, path):
    """ansys_rst_cyclic: the base sector repeated round the cyclic axis (a local
    coordinate system's Z in cyc12.rst), as pymapdl-reader's full rotor."""
    ref = np.load(RST / "pymapdl_reference.npz")
    key = path.name
    rotor = ref[f"{key}/rotor_points"]
    base = py_rst._Model(py_rst._open(str(path)), True)
    r = base.files[0]
    if r.cs_cord > 1:
        axes, origin = r.coordinate_system(r.cs_cord)
        axis = axes[2] / np.linalg.norm(axes[2])
    else:
        axis, origin = np.array([0.0, 0.0, 1.0]), np.zeros(3)
    # The rotor holds the base sector's points (those of its cells, in point
    # order) once per sector, sector by sector.
    used = set()
    for e, loc in zip(base.deck["elements"], base.locs):
        if loc is not None and e["id"] <= r.cs_els:
            used.update(int(v) for v in base.mesh.cells[loc[0]].data[loc[1]])
    rank = {p: k for k, p in enumerate(sorted(used))}
    for step in range(len(ref[f"{key}/times"])):
        mesh = engine.read_cyclic(path, time_step=step, lenient=True)
        assert mesh.field_data["ansys:sectors"].tolist() == [r.n_sectors]
        assert len(mesh.points) == len(rotor) == r.n_sectors * len(used)
        scale = np.abs(rotor).max()
        np.testing.assert_allclose(
            np.sort(mesh.points.round(9), axis=0),
            np.sort(rotor.round(9), axis=0),
            atol=1e-8 * scale,
        )
        numbers = ref[f"{key}/{step}/rotor_nnum"]
        values = ref[f"{key}/{step}/rotor_values"]
        stress = ref[f"{key}/{step}/rotor_stress"]
        base_rows = [rank[base.node_index[int(n)]] for n in numbers]
        xyz = base.mesh.points[[base.node_index[int(n)] for n in numbers]] - origin
        for i in range(r.n_sectors):
            rows = [i * len(used) + k for k in base_rows]
            q = np.array(py_rst._axis_rotation(axis, 2 * math.pi * i / r.n_sectors))
            np.testing.assert_allclose(
                mesh.points[rows], xyz @ q.reshape(3, 3).T + origin, atol=1e-12 * scale
            )
            np.testing.assert_allclose(
                mesh.point_data["U"][rows],
                values[i][:, :3],
                rtol=0,
                atol=1e-12 * np.abs(values).max(),
            )
            np.testing.assert_allclose(
                mesh.point_data["S"][rows],
                stress[i],
                rtol=0,
                atol=1e-12 * np.nanmax(np.abs(stress)),
            )
        sector = np.concatenate(mesh.cell_data["ansys:sector"])
        assert sorted(set(sector.tolist())) == list(range(r.n_sectors))


def test_cyclic_refusals_and_dispatch(engine, tmp_path):
    with pytest.raises(meshioplusplus.ReadError, match="not a cyclic-symmetry model"):
        engine.read_cyclic(RST / "file.rst")
    # A cyclic model of another analysis type (here: transient, 4) is refused.
    other = tmp_path / "transient.rst"
    modal_variant(RST / "cyc12.rst", other, (0, 0, 0), (1.0, 2.0, 3.0))
    data = bytearray(other.read_bytes())
    header = int(np.frombuffer(bytes(data[:4]), "<i4")[0]) + 3
    data[(header + 9) * 4 : (header + 10) * 4] = np.int32(4).tobytes()
    other.write_bytes(bytes(data))
    with pytest.raises(meshioplusplus.ReadError, match="only static and modal"):
        engine.read_cyclic(other)
    mesh = meshioplusplus.read(RST / "cyclic_v182.rst", file_format="ansys_rst_cyclic")
    assert mesh.field_data["ansys:sectors"].tolist() == [15]
    assert meshioplusplus.read_metadata(
        RST / "cyc12.rst", file_format="ansys_rst_cyclic"
    )["time_values"] == pytest.approx(py_rst.time_values(RST / "cyc12.rst"))


def modal_variant(source, target, harmonic, times):
    """Write ``source`` (a static cyclic results file) to ``target`` as a modal
    one: analysis type 2, the harmonic index table and the set times rewritten.
    ``tools/gen_ansys_rst_reference.py`` holds the same function."""
    data = bytearray(pathlib.Path(source).read_bytes())
    words = np.frombuffer(bytes(data), "<i4")
    header = int(words[0]) + 3  # the result header follows the standard one
    h = words[header + 2 : header + 2 + int(words[header])]

    def pointer(lo, hi):
        return (int(h[lo]) & 0xFFFFFFFF) | (int(h[hi]) << 32)

    def put(word, value):
        raw = value.tobytes()
        data[word * 4 : word * 4 + len(raw)] = raw

    put(header + 2 + 7, np.int32(2))
    cyc = pointer(16, 43)
    for k, v in enumerate(harmonic):
        put(cyc + 2 + k, np.int32(v))
    tim = pointer(11, 41)
    for k, v in enumerate(times):
        put(tim + 2 + 2 * k, np.float64(v))
    pathlib.Path(target).write_bytes(bytes(data))


def _rotor_rows(path):
    """The base model, its cyclic axis and origin, and each base point's row in
    a sector of the full rotor (the base sector's points in point order)."""
    base = py_rst._Model(py_rst._open(str(path)), True)
    r = base.files[0]
    if r.cs_cord > 1:
        axes, origin = r.coordinate_system(r.cs_cord)
        axis = axes[2] / np.linalg.norm(axes[2])
    else:
        axis, origin = np.array([0.0, 0.0, 1.0]), np.zeros(3)
    used = set()
    for e, loc in zip(base.deck["elements"], base.locs):
        if loc is not None and e["id"] <= r.cs_els:
            used.update(int(v) for v in base.mesh.cells[loc[0]].data[loc[1]])
    rank = {p: k for k, p in enumerate(sorted(used))}
    return base, axis, origin, rank


def _sector_rotation(axis, n, i):
    return np.array(py_rst._axis_rotation(axis, 2 * math.pi * i / n)).reshape(3, 3)


def test_modal_cyclic_mode_pair_matches_pymapdl_reader(engine, tmp_path):
    """A modal cyclic set is combined with its mode pair (the neighbouring set
    of the same frequency) per sector, as pymapdl-reader's CyclicResult does:
    the stresses match its full rotor on every sector. pymapdl-reader rotates
    modal displacements about the global Z axis whatever the cyclic axis
    (cyc12.rst's is a local system's), so they are compared unrotated."""
    path = tmp_path / "cyc12_modal_pair.rst"
    modal_variant(RST / "cyc12.rst", path, (4, -4, 9), (5.0, 5.0, 7.0))
    ref = np.load(RST / "pymapdl_reference.npz")
    key = path.name
    base, axis, _, rank = _rotor_rows(path)
    n = base.files[0].n_sectors
    for step, harmonic in enumerate((4, -4, 9)):
        mesh = engine.read_cyclic(path, time_step=step, lenient=True)
        assert mesh.field_data["ansys:harmonic_index"].tolist() == [harmonic]
        numbers = ref[f"{key}/{step}/rotor_nnum"]
        values = ref[f"{key}/{step}/rotor_values"]
        stress = ref[f"{key}/{step}/rotor_stress"]
        rows = [rank[base.node_index[int(v)]] for v in numbers]
        for i in range(n):
            here = [i * len(rank) + k for k in rows]
            ours = mesh.point_data["U"][here] @ _sector_rotation(axis, n, i)
            theirs = values[i][:, :3] @ _sector_rotation([0.0, 0.0, 1.0], n, i)
            np.testing.assert_allclose(
                ours, theirs, rtol=0, atol=1e-12 * np.abs(values).max()
            )
            np.testing.assert_allclose(
                mesh.point_data["S"][here],
                stress[i],
                rtol=0,
                atol=1e-12 * np.nanmax(np.abs(stress)),
            )


def test_modal_cyclic_duplicate_sector(engine, tmp_path):
    """Without a mode pair, a modal set's other half is its duplicate sector
    (nodes past csNds, paired with the base sector's in number order): sector i
    is scale * (x cos(h theta_i) - x' sin(h theta_i)), rotated. The base sector
    matches pymapdl-reader (whose own duplicate-sector branch never runs: it
    tests the reduced array's size against the full node count)."""
    path = tmp_path / "cyc12_modal_dup.rst"
    modal_variant(RST / "cyc12.rst", path, (2, 3, 0), (1.0 / 3.0, 2.0 / 3.0, 1.0))
    ref = np.load(RST / "pymapdl_reference.npz")
    key = path.name
    base, axis, _, rank = _rotor_rows(path)
    r = base.files[0]
    n = r.n_sectors
    pairs = sorted(v for v in base.node_index if v <= r.cs_nds)
    dups = sorted(v for v in base.node_index if v > r.cs_nds)
    assert len(pairs) == len(dups)
    for step, harmonic in enumerate((2, 3, 0)):
        plain = engine.read(path, time_step=step, lenient=True)
        mesh = engine.read_cyclic(path, time_step=step, lenient=True)
        u = plain.point_data["U"]
        x = u[[base.node_index[v] for v in pairs]]
        x_dup = u[[base.node_index[v] for v in dups]]
        rows = [rank[base.node_index[v]] for v in pairs]
        single = harmonic == 0 or 2 * abs(harmonic) == n
        scale = 1 / math.sqrt(n) if single else 1 / math.sqrt(n / 2)
        for i in range(n):
            phase = 2 * math.pi * harmonic * i / n
            expected = scale * (x * math.cos(phase) - x_dup * math.sin(phase))
            here = [i * len(rank) + k for k in rows]
            np.testing.assert_allclose(
                mesh.point_data["U"][here] @ _sector_rotation(axis, n, i),
                expected,
                rtol=0,
                atol=1e-14 * np.abs(u).max(),
            )
        numbers = ref[f"{key}/{step}/rotor_nnum"]
        rows = [rank[base.node_index[int(v)]] for v in numbers]
        np.testing.assert_allclose(
            mesh.point_data["U"][rows], ref[f"{key}/{step}/rotor_values"][0][:, :3]
        )
        np.testing.assert_allclose(
            mesh.point_data["S"][rows], ref[f"{key}/{step}/rotor_stress"][0]
        )


_RAW = sorted(
    {
        (k.split("/raw/")[0], k.split("/raw/")[1].split("/")[0])
        for k in np.load(RST / "pymapdl_reference.npz").files
        if "/raw/" in k
    }
)


@pytest.mark.parametrize("key,name", _RAW, ids=[f"{k}:{n}" for k, n in _RAW])
def test_raw_element_records_match_pymapdl_reader(engine, key, name):
    """ENG, EMS, EMN, EGR, EFX, ENL, EPT ...: one row per element, the record as
    written (NaN-padded to the longest), as pymapdl-reader's
    ``element_solution_data`` reads it. On the release 13 files pymapdl-reader
    reads twice each record's length (it counts 4-byte words as values): only
    the first half is the record."""
    path = RST / key
    ref = np.load(RST / "pymapdl_reference.npz")
    enum = ref[f"{key}/raw/{name}/enum"]
    lengths = ref[f"{key}/raw/{name}/lengths"]
    values = np.split(ref[f"{key}/raw/{name}/values"], np.cumsum(lengths)[:-1])
    if key in ("beam44.rst", "temp_v13.rst"):
        lengths = lengths // 2
    model = py_rst._Model(py_rst._open(str(path)), True)
    where = {
        e["id"]: loc
        for e, loc in zip(model.deck["elements"], model.locs)
        if loc is not None
    }
    mesh = engine.read(path, lenient=True)
    data = mesh.cell_data[name]
    checked = 0
    for number, length, expected in zip(enum, lengths, values):
        if int(number) not in where or not length:
            continue
        b, row = where[int(number)][:2]
        got = data[b][row]
        np.testing.assert_array_equal(got[:length], expected[:length])
        assert np.isnan(got[length:]).all()
        checked += 1
    assert checked


def test_zero_records(engine):
    """A negative pointer-table entry ``-n`` stands for ``n`` zeros: MAPDL does
    not write an all-zero record (cyclic_v182.rst's unloaded bricks)."""
    mesh = engine.read(RST / "cyclic_v182.rst")
    cells = np.concatenate([c.reshape(-1, 6) for c in mesh.cell_data["S"]])
    corners = cells[np.isfinite(cells).all(1)]
    assert (corners == 0).all(1).sum() >= 8 and (corners != 0).any()


def test_modal_file(engine):
    """file.rst: a modal analysis, six modes of 40 SOLID186 bricks."""
    path = RST / "file.rst"
    np.testing.assert_allclose(
        engine.time_values(path),
        [7366.495, 7366.495, 11504.895, 17285.705, 17285.705, 20137.193],
        atol=1e-3,
    )
    mesh = engine.read(path, time_step=-1)
    assert [b.type for b in mesh.cells] == ["hexahedron20"]
    assert mesh.cells[0].data.shape == (40, 20)
    assert mesh.field_data["meshio:time"][0] == pytest.approx(20137.193, abs=1e-3)
    assert mesh.field_data["ansys:cumulative"].tolist() == [6]
    first = engine.read(path)
    np.testing.assert_allclose(
        first.point_data["U"][0], [28.962, -28.248, -0.309], atol=1e-3
    )
    np.testing.assert_array_equal(first.cell_data["ansys:element"][0], 186)


def test_thermal_file_has_temperatures(engine):
    mesh = engine.read(RST / "file.rth")
    assert sorted(mesh.point_data) == ["RF_TEMP", "TEMP"]
    assert np.isfinite(mesh.point_data["TEMP"]).all()


def test_components_are_regions(engine):
    mesh = engine.read(RST / "hex_201.rst")
    names = {(r.kind, r.name): len(r.entries) for r in mesh.regions}
    assert names[("point", "NCOMP2")] == 98
    assert names[("cell", "ECOMP1")] == 22


def test_cyclic_model_reads_its_sector_with_a_warning(engine, capfd):
    mesh = engine.read(RST / "cyclic_v182.rst")
    assert "cyclic-symmetry" in capfd.readouterr().err
    assert len(mesh.points) == 460


def test_selective_reads(engine):
    path = RST / "shell181_2021R1.rst"
    assert not engine.read(path, points_only=True).point_data
    assert sorted(engine.read(path, arrays=["ROT"]).point_data) == ["ROT"]
    assert not engine.read(path, arrays=[]).point_data
    only = engine.read(path, arrays=["S", "RF"])
    assert sorted(only.point_data) == ["RF", "S"]
    assert sorted(k for k in only.cell_data if not k.startswith("ansys:")) == ["S"]
    with pytest.raises(meshioplusplus.ReadError, match="out of range"):
        engine.read(path, time_step=1)


def test_sequences_and_metadata(tmp_path):
    path = RST / "file.rst"
    assert len(meshioplusplus.read_metadata(path)["time_values"]) == 6
    meshioplusplus.write_sequence(
        str(tmp_path / "mode_{step}.vtu"), meshioplusplus.read_sequence(path)
    )
    written = sorted(tmp_path.glob("mode_*.vtu"))
    assert len(written) == 6
    np.testing.assert_allclose(
        meshioplusplus.read(written[-1]).point_data["U"],
        meshioplusplus.read(path, time_step=5).point_data["U"],
    )


def test_registered_and_sniffed(tmp_path):
    out = meshioplusplus.formats()
    assert out["extensions"][".rst"] == ["ansys_rst"]
    assert out["extensions"][".rth"] == ["ansys_rst"]
    assert "ansys_rst" in out["readable"] and "ansys_rst" not in out["writable"]
    assert "ansys_rst_cyclic" in out["readable"]
    assert "ansys_rst_cyclic" not in out["writable"]
    renamed = tmp_path / "results.bin"
    renamed.write_bytes((RST / "beam44.rst").read_bytes())
    assert meshioplusplus.sniff_format(renamed) == "ansys_rst"
    assert hasattr(_core, "ansys_rst_read")


# -- a synthetic results file ------------------------------------------------------


def _i32(value):
    return value - (1 << 32) if value >= 1 << 31 else value


class _Writer:
    """Records as MAPDL lays them out: length in words, flags, data, trailer."""

    def __init__(self):
        self.words = []

    def raw(self, data_words, flags):
        ptr = len(self.words)
        self.words += [len(data_words), _i32(flags << 24), *data_words, 0]
        return ptr

    def ints(self, values, size=None):
        values = list(values) + [0] * ((size or len(values)) - len(values))
        return self.raw(values, 0x80)

    def doubles(self, values):
        return self.raw(np.asarray(values, "<f8").view("<i4").tolist(), 0x00)

    def shorts(self, values):
        values = list(values) + [0] * (len(values) % 2)
        return self.raw(np.asarray(values, "<i2").view("<i4").tolist(), 0xC0)

    def bsparse_doubles(self, values):
        bits = sum(1 << k for k, v in enumerate(values) if v != 0)
        packed = np.asarray([v for v in values if v != 0], "<f8").view("<i4").tolist()
        return self.raw([len(values), _i32(bits)] + packed, 0x08)

    def patch(self, ptr, k, value):
        self.words[ptr + 2 + k] = _i32(value & 0xFFFFFFFF)

    def save(self, path):
        pathlib.Path(path).write_bytes(np.asarray(self.words, "<i4").tobytes())


_ANGLES = (30.0, 45.0, 60.0)  # node 3: THXY, THYZ, THZX
_NEQV = [8, 7, 6, 5, 4, 3, 2, 1]  # the solution's node order
_HEX = [
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
]


# The synthetic element solution: the brick's stress at its node i (element
# axes, turned 90 degrees about Z to the global ones), the shell's bottom and
# top stresses, and three reactions (node 3's UX in its rotated axes).
def _brick_stress(i):
    return [1.0 + i, 2.0, 3.0, 0.5, 0.25, 0.125]


def _shell_stress(i, top):
    return [10.0 + i + 100.0 * top, 20.0, 0.0, 1.0, 0.0, 0.0]


_REACTIONS = [(3, 1, 7.0), (1, 3, -2.0), (2, 4, 5.0)]  # node, DOF position, value


def _write_synthetic(path, zlib=False, global_nnod=0, sectors=1, elements=False):
    w = _Writer()
    w.ints([12], 100)  # the standard header: file 12
    header = w.ints([12, 8, 8, 10, 4, 1, 1, 0, 2], 80)
    w.patch(header, 20, sectors)
    w.patch(header, 48, global_nnod)
    w.patch(header, 39, 1)  # rstsprs: ENS holds six items per node
    w.patch(header, 14, w.ints(_NEQV))

    geometry = w.ints([0, 2 if elements else 1, 0, 8, 2 if elements else 1], 80)
    w.patch(header, 15, geometry)
    ety = w.ints([0, 0] if elements else [0])
    w.patch(geometry, 20, ety)
    solid185 = w.ints([1, 185, 0], 100)
    w.patch(solid185, 60, 8)  # nodelm
    w.patch(solid185, 62, 8)  # nodfor
    w.patch(solid185, 93, 8)  # nodstr
    w.patch(ety, 0, solid185 - ety)
    if elements:
        shell181 = w.ints([2, 181, 0], 100)  # KEYOPT(8) = 0: bottom and top
        w.patch(shell181, 60, 4)
        w.patch(shell181, 62, 4)
        w.patch(shell181, 93, 4)
        w.patch(ety, 1, shell181 - ety)
    nodes = None
    for number, xyz in enumerate(_HEX, start=1):
        angles = _ANGLES if number == 3 else (0.0, 0.0, 0.0)
        values = [float(number), *map(float, xyz), *angles]
        ptr = w.bsparse_doubles(values) if number == 1 else w.doubles(values)
        nodes = ptr if nodes is None else nodes
    w.patch(geometry, 26, nodes)
    eid = w.ints([0, 0, 0, 0] if elements else [0, 0])
    w.patch(geometry, 28, eid)
    w.patch(eid, 0, w.shorts([2, 1, 3, 4, 0, 1, 0, 0, 11, 0, *range(1, 9)]) - eid)
    if elements:
        w.patch(eid, 2, w.ints([2, 2, 3, 4, 0, 1, 0, 0, 12, 0, 1, 2, 3, 4]) - eid)
    comp = w.ints([1] + [int.from_bytes(b"FACE", "big")] + [0x20202020] * 7 + [1, -4])
    w.patch(geometry, 50, comp)
    w.patch(geometry, 48, 1)

    sets = []
    for s in range(2):
        base = w.ints([0, 1, 8], 150)
        w.patch(base, 19, 4)
        for k, code in enumerate([1, 2, 3, 20]):
            w.patch(base, 20 + k, code)
        if s == 0:  # every node: U = (n, 0, 0), TEMP = 1000 + n
            rows = [[n, 0.0, 0.0, 1000.0 + n] for n in _NEQV]
            nsl = w.doubles(np.ravel(rows))
        else:  # nodes 5 and 2 only, node 5's TEMP undefined
            rows = [[0.0, 5.0, 0.0, 2.0**100], [0.0, 2.0, 0.0, 1002.0]]
            nsl = w.doubles(np.ravel(rows))
            w.ints([_NEQV.index(5) + 1, _NEQV.index(2) + 1])
        w.patch(base, 104, nsl - base)
        if elements and s == 0:
            _write_element_solution(w, base)
        sets.append(base)
    dsi = w.ints(sets, 20)
    tim = w.doubles([0.5, 1.0] + [0.0] * 8)
    lsp = w.ints([1, 1, 1, 1, 2, 2], 30)
    for k, ptr in ((10, dsi), (11, tim), (12, lsp)):
        w.patch(header, k, ptr)
    if zlib:
        w.words[tim + 1] = _i32(0x20 << 24)
    w.save(path)


def _write_element_solution(w, base):
    """Set 0's reactions and element solution (``ptrRF``, ``ptrESL``)."""
    table = []
    for node, k, _ in _REACTIONS:
        index = _NEQV.index(node) * 4 + k  # (N - 1) * numdof + k, N 1-based
        table += [index, 0]
    rf = w.ints(table)
    w.doubles([v for _, _, v in _REACTIONS])
    w.patch(base, 7, len(_REACTIONS))
    w.patch(base, 106, rf - base)

    esl = w.ints([0, 0, 0, 0])
    w.patch(base, 118, esl - base)
    # The brick: stresses, a 90-degree element rotation, no elastic strain (an
    # all-zero record, entry -56) and nodal forces.
    brick = w.ints([0] * 25)
    w.patch(esl, 0, brick - esl)
    w.patch(brick, 2, w.doubles(np.ravel([_brick_stress(i) for i in range(8)])) - brick)
    w.patch(brick, 5, -56)
    w.patch(brick, 9, w.doubles([90.0, 0.0, 0.0]) - brick)
    forces = [[i, 0.0, 0.0, -i] for i in range(8)]
    w.patch(brick, 1, w.doubles(np.ravel(forces)) - brick)
    # Energies (ENG), read as written: the undefined value becomes NaN.
    w.patch(brick, 3, w.doubles([1.0, 2.0, 2.0**100, 4.0]) - brick)
    # The shell: bottom then top, four corners each.
    shell = w.ints([0] * 25)
    w.patch(esl, 2, shell - esl)
    w.patch(shell, 3, -2)  # an all-zero energy record, not written
    stresses = [_shell_stress(i, 0) for i in range(4)] + [
        _shell_stress(i, 1) for i in range(4)
    ]
    w.patch(shell, 2, w.doubles(np.ravel(stresses)) - shell)


def _rotation(xy, yz, zx):
    def about(axis, degrees):
        c, s = math.cos(math.radians(degrees)), math.sin(math.radians(degrees))
        i, j = [(1, 2), (2, 0), (0, 1)][axis]
        m = np.eye(3)
        m[i, i], m[i, j], m[j, i], m[j, j] = c, -s, s, c
        return m

    return about(2, xy) @ about(0, yz) @ about(1, zx)


def test_synthetic_file(engine, tmp_path):
    path = tmp_path / "synthetic.rst"
    _write_synthetic(path)
    assert engine.time_values(path) == [0.5, 1.0]

    mesh = engine.read(path)
    assert [b.type for b in mesh.cells] == ["hexahedron"]
    np.testing.assert_array_equal(mesh.cells[0].data, [list(range(8))])
    np.testing.assert_array_equal(mesh.points, np.array(_HEX, dtype=float))
    for name, value in (("ansys:mat", 2), ("ansys:real", 3), ("ansys:secnum", 4)):
        assert mesh.cell_data[name][0].tolist() == [value]
    assert [(r.kind, r.name, r.entries.tolist()) for r in mesh.regions] == [
        ("point", "FACE", [0, 1, 2, 3])
    ]
    assert mesh.field_data["meshio:time"].tolist() == [0.5]
    np.testing.assert_array_equal(mesh.point_data["TEMP"], 1000.0 + np.arange(1, 9))
    expected = np.array([[n, 0.0, 0.0] for n in range(1, 9)])
    # Node 3's axes are turned about Z, then the new X, then the newest Y.
    expected[2] = _rotation(*_ANGLES) @ expected[2]
    np.testing.assert_allclose(mesh.point_data["U"], expected, rtol=1e-14, atol=1e-14)

    partial = engine.read(path, time_step=1)
    assert partial.field_data["ansys:substep"].tolist() == [2]
    temp, u = partial.point_data["TEMP"], partial.point_data["U"]
    assert np.isnan(temp[[0, 2, 3, 4, 5, 6, 7]]).all()  # node 5's is undefined
    assert temp[1] == 1002.0
    np.testing.assert_array_equal(u[1], [0.0, 2.0, 0.0])
    np.testing.assert_array_equal(u[4], [0.0, 5.0, 0.0])
    assert np.isnan(u[[0, 2, 3, 5, 6, 7]]).all()


def test_synthetic_element_solution(engine, tmp_path):
    """A rotated brick, a layered shell, an all-zero record and reactions."""
    path = tmp_path / "elements.rst"
    _write_synthetic(path, elements=True)
    mesh = engine.read(path)
    assert [b.type for b in mesh.cells] == ["hexahedron", "quad"]

    # The brick's element x axis is the global y axis.
    q = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])

    def to_global(t):
        xx, yy, zz, xy, yz, xz = t
        m = q @ np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]) @ q.T
        return [m[0, 0], m[1, 1], m[2, 2], m[0, 1], m[1, 2], m[0, 2]]

    # Every block's arrays have the brick's eight nodes; the quad's last four
    # are NaN.
    np.testing.assert_array_equal(mesh.field_data["ansys:layout:S"], [8, 6])

    def cell(name, block):
        return mesh.cell_data[name][block][0].reshape(8, -1)

    brick = np.array([to_global(_brick_stress(i)) for i in range(8)])
    np.testing.assert_allclose(cell("S", 0), brick, atol=1e-12)
    np.testing.assert_array_equal(cell("EPEL", 0), np.zeros((8, 6)))
    assert np.isnan(mesh.cell_data["EPEL"][1]).all()  # the shell has none
    shell_bottom = np.array([_shell_stress(i, 0) for i in range(4)])
    shell_top = np.array([_shell_stress(i, 1) for i in range(4)])
    np.testing.assert_array_equal(cell("S", 1)[:4], shell_bottom)
    assert np.isnan(cell("S", 1)[4:]).all()
    np.testing.assert_array_equal(cell("S@top", 1)[:4], shell_top)
    assert np.isnan(mesh.cell_data["S@top"][0]).all()

    # Averaged at the nodes: the bottom face (nodes 1-4) with the shell's bottom.
    expected = brick.copy()
    expected[:4] = (brick[:4] + shell_bottom) / 2
    np.testing.assert_allclose(mesh.point_data["S"], expected, atol=1e-12)
    top = np.full((8, 6), np.nan)
    top[:4] = shell_top
    np.testing.assert_array_equal(mesh.point_data["S@top"], top)
    np.testing.assert_array_equal(mesh.point_data["EPEL"][4:], np.zeros((4, 6)))

    np.testing.assert_array_equal(mesh.field_data["ansys:layout:ENF"], [8, 4])
    forces = cell("ENF", 0)
    np.testing.assert_array_equal(forces, [[i, 0.0, 0.0, -i] for i in range(8)])
    assert np.isnan(mesh.cell_data["ENF"][1]).all()
    # Raw records: one row per element, NaN-padded to the longest.
    np.testing.assert_array_equal(mesh.cell_data["ENG"][0], [[1.0, 2.0, np.nan, 4.0]])
    np.testing.assert_array_equal(
        mesh.cell_data["ENG"][1], [[0.0, 0.0, np.nan, np.nan]]
    )
    assert "EMS" not in mesh.cell_data

    rf = mesh.point_data["RF"]
    np.testing.assert_allclose(rf[2], _rotation(*_ANGLES) @ [7.0, 0.0, 0.0], atol=1e-14)
    np.testing.assert_array_equal(rf[0], [0.0, 0.0, -2.0])
    assert np.isnan(rf[[1, 3, 4, 5, 6, 7]]).all()
    temp = mesh.point_data["RF_TEMP"]
    assert temp[1] == 5.0 and np.isnan(np.delete(temp, 1)).all()

    later = engine.read(path, time_step=1)  # set 1 has no element solution
    assert "S" not in later.point_data and "RF" not in later.point_data


@pytest.mark.parametrize(
    "options, match",
    [
        ({"zlib": True}, "zlib"),
        ({"global_nnod": 16}, "distributed"),
    ],
)
def test_refusals(engine, tmp_path, options, match):
    path = tmp_path / "bad.rst"
    _write_synthetic(path, **options)
    with pytest.raises(meshioplusplus.ReadError, match=match):
        engine.read(path)
    junk = tmp_path / "junk.rst"
    junk.write_bytes(b"\x00" * 64)
    with pytest.raises(meshioplusplus.ReadError, match="not a MAPDL results file"):
        engine.read(junk)
    cut = tmp_path / "cut.rst"
    cut.write_bytes((RST / "beam44.rst").read_bytes()[:5000])
    with pytest.raises(meshioplusplus.ReadError, match="outside the file|past the end"):
        engine.read(cut)


@pytest.mark.parametrize("ext", ["vtu", "xdmf"])
def test_element_results_are_writable(tmp_path, ext):
    """Per-element-node arrays are rectangular and as wide in every block
    (quad8 and hexahedron20 here), so the writers hold them: they round-trip."""
    mesh = meshioplusplus.read(DIST / "file0.rst")
    assert len({a.shape[1] for a in mesh.cell_data["S"]}) == 1
    target = tmp_path / f"out.{ext}"
    meshioplusplus.write(target, mesh)
    back = meshioplusplus.read(target)
    for name in ("S", "EPEL", "ENF"):
        for a, b in zip(mesh.cell_data[name], back.cell_data[name]):
            np.testing.assert_array_equal(a, b)
