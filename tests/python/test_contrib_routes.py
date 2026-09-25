"""The vendor-runtime routes under ``contrib/``: each script runs here against a
stand-in of its vendor API (only the calls it makes, as the vendor documents
them), and its output is read back through meshio++.

These prove the scripts' own logic -- element mapping, sets, field positions,
times, the files written -- not the vendor APIs themselves: the run against
each licensed tool is listed in doc/roadmap.md ("Awaiting a licensed run")."""

import importlib.util
import pathlib
import sys
import types
from collections import OrderedDict

import numpy as np
import pytest

import meshioplusplus

CONTRIB = pathlib.Path(__file__).resolve().parents[2] / "contrib"


def _load(relative, name):
    spec = importlib.util.spec_from_file_location(name, CONTRIB / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _ns(**kwargs):
    return types.SimpleNamespace(**kwargs)


# --------------------------------------------------------------------------- #
# The shared series writer                                                    #
# --------------------------------------------------------------------------- #


def test_series_writer_round_trips(tmp_path):
    series_mod = _load("common/mio_vtu_series.py", "mio_vtu_series_test")
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], float)
    with series_mod.PvdSeries(str(tmp_path / "run.pvd")) as series:
        for t in (0.0, 0.5):
            series.add(
                t,
                points,
                [(10, [[0, 1, 2, 3]]), (5, [[1, 2, 4]])],
                point_data={"T": points[:, 0] + t},
                cell_data={"id": np.array([7, 8], dtype=np.int32)},
                part=0,
                name="body",
            )
    steps = list(meshioplusplus.read_sequence(tmp_path / "run.pvd"))
    assert [t for t, _ in steps] == [0.0, 0.5]
    mesh = steps[1][1]
    assert [c.type for c in mesh.cells] == ["tetra", "triangle"]
    np.testing.assert_array_equal(mesh.point_data["T"], points[:, 0] + 0.5)
    assert [list(a) for a in mesh.cell_data["id"]] == [[7], [8]]


# --------------------------------------------------------------------------- #
# Abaqus .odb                                                                 #
# --------------------------------------------------------------------------- #


class _Repo(OrderedDict):
    """An ODB repository: a mapping with keys() in creation order."""


def _abaqus_constants():
    return _ns(
        NODAL="NODAL",
        CENTROID="CENTROID",
        ELEMENT_NODAL="ELEMENT_NODAL",
        INTEGRATION_POINT="INTEGRATION_POINT",
        WHOLE_ELEMENT="WHOLE_ELEMENT",
    )


class _Field:
    def __init__(self, name, position, values):
        self.name = name
        self.locations = [_ns(position=position)]
        self._values = values  # position -> list of value objects

    def getSubset(self, region, position):
        return _ns(values=self._values[position])


def _odb():
    """Two instances: a C3D4 tetrahedron (+ a MASS element with no cell) and an
    S4R shell; a static step with two frames and a frequency step."""
    tet_nodes = [
        _ns(label=10 + i, coordinates=xyz)
        for i, xyz in enumerate([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
    ]
    tet = _ns(
        name="SOLID-1",
        nodes=tet_nodes,
        elements=[
            _ns(label=1, type="C3D4", connectivity=(10, 11, 12, 13)),
            _ns(label=2, type="MASS", connectivity=(10,)),
        ],
        nodeSets=_Repo(BASE=_ns(nodes=[tet_nodes[0], tet_nodes[1]])),
        elementSets=_Repo(),
    )
    shell_nodes = [
        _ns(label=1 + i, coordinates=xyz)
        for i, xyz in enumerate([(2, 0, 0), (3, 0, 0), (3, 1, 0), (2, 1, 0)])
    ]
    shell = _ns(
        name="SKIN-1",
        nodes=shell_nodes,
        elements=[_ns(label=5, type="S4R", connectivity=(1, 2, 3, 4))],
        nodeSets=_Repo(),
        elementSets=_Repo(ALL=_ns(elements=[_ns(label=5)])),
    )
    assembly = _ns(
        instances=_Repo([("SOLID-1", tet), ("SKIN-1", shell)]),
        nodeSets=_Repo(
            TIP=_ns(
                instanceNames=("SOLID-1", "SKIN-1"),
                nodes=[[tet_nodes[3]], [shell_nodes[2]]],
            )
        ),
        elementSets=_Repo(),
    )

    def frame(value, scale, domain="TIME"):
        u_values = [
            _ns(nodeLabel=n.label, data=(scale * n.label, 0.0, 0.0))
            for n in tet_nodes + shell_nodes
        ]
        s_values = [
            _ns(elementLabel=1, data=(scale, 0, 0, 0, 0, 0)),
            _ns(elementLabel=1, data=(3 * scale, 0, 0, 0, 0, 0)),
            _ns(elementLabel=5, data=(2 * scale, 0, 0, 0, 0, 0)),
        ]
        return _ns(
            frameValue=value,
            domain=domain,
            fieldOutputs=_Repo(
                U=_Field("U", "NODAL", {"NODAL": u_values}),
                S=_Field("S", "INTEGRATION_POINT", {"CENTROID": s_values}),
            ),
        )

    steps = _Repo(
        [
            ("Load", _ns(totalTime=0.0, frames=[frame(0.0, 0.0), frame(1.0, 1.0)])),
            (
                "Modes",
                _ns(totalTime=1.0, frames=[frame(0.0, 5.0, "FREQUENCY")]),
            ),
        ]
    )
    return _ns(rootAssembly=assembly, steps=steps)


def test_abaqus_odb_route(tmp_path, capsys):
    script = _load("abaqus_odb/odb_to_vtu.py", "odb_to_vtu_test")
    out = tmp_path / "job.pvd"
    written = script.export(_odb(), str(out), constants=_abaqus_constants())
    assert written == 3
    assert "MASS" in capsys.readouterr().out  # the element with no cell

    steps = list(meshioplusplus.read_sequence(out))
    assert [t for t, _ in steps] == [0.0, 1.0, 2.0]  # the modal frame after 1.0
    mesh = steps[1][1]
    assert sorted(r.name for r in mesh.regions) == ["SKIN-1", "SOLID-1"]
    assert [c.type for c in mesh.cells] == ["tetra", "quad"]
    # U by node label, S averaged over the element's integration points.
    np.testing.assert_array_equal(
        mesh.point_data["U"][:, 0], [10, 11, 12, 13, 1, 2, 3, 4]
    )
    assert [a[0, 0] for a in mesh.cell_data["S"]] == [2.0, 2.0]
    np.testing.assert_array_equal(mesh.point_data["set:BASE"], [1, 1, 0, 0, 0, 0, 0, 0])
    np.testing.assert_array_equal(mesh.point_data["set:TIP"], [0, 0, 0, 1, 0, 0, 1, 0])
    assert [list(a) for a in mesh.cell_data["set:ALL"]] == [[0], [1]]
    assert [list(a) for a in mesh.cell_data["abaqus:element_label"]] == [[1], [5]]


def test_abaqus_element_table():
    script = _load("abaqus_odb/odb_to_vtu.py", "odb_to_vtu_table")
    assert script.vtk_type("C3D20R", 20) == 25
    assert script.vtk_type("C3D10M", 10) == 24
    assert script.vtk_type("CPS8R", 8) == 23
    assert script.vtk_type("S8R5", 8) == 23
    assert script.vtk_type("SC8R", 8) == 12  # a continuum shell is a solid
    assert script.vtk_type("B32", 3) == 21
    assert script.vtk_type("SPRING1", 1) is None


# --------------------------------------------------------------------------- #
# Marc .t16 (PyPost)                                                          #
# --------------------------------------------------------------------------- #


class _PyPost:
    """Position 0 is the model; positions 1 and 2 are increments at t = 0.5, 1."""

    def __init__(self):
        self.position = 0
        self.xyz = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (2, 0, 0)]
        self.els = [
            _ns(type=11, items=[1, 2, 3, 4]),  # quad
            _ns(type=6, items=[2, 5, 3]),  # triangle
            _ns(type=9999, items=[1, 5]),  # unknown
        ]
        self.time = 0.0

    def increments(self):
        return 3

    def moveto(self, position):
        self.position = position
        self.time = 0.5 * position

    def nodes(self):
        return len(self.xyz)

    def node(self, i):
        x, y, z = self.xyz[i]
        return _ns(x=x, y=y, z=z)

    def node_id(self, i):
        return i + 1

    def elements(self):
        return len(self.els)

    def element(self, e):
        return self.els[e]

    def element_id(self, e):
        return 100 + e

    def node_scalars(self):
        return 1

    def node_scalar_label(self, k):
        return "Temperature"

    def node_scalar(self, i, k):
        return 10.0 * self.position + i

    def node_vectors(self):
        return 1

    def node_vector_label(self, k):
        return "Displacement"

    def node_vector(self, i, k):
        return _ns(x=self.position * i, y=0.0, z=0.0)

    def element_scalars(self):
        return 1

    def element_scalar_label(self, k):
        return "Equivalent Von Mises Stress"

    def element_scalar(self, e, k):
        return [_ns(id=n, value=float(n)) for n in self.els[e].items]

    def element_tensors(self):
        return 1

    def element_tensor_label(self, k):
        return "Stress"

    def element_tensor(self, e, k):
        return [
            _ns(t11=1.0 * e, t22=2.0, t33=3.0, t12=4.0, t23=5.0, t13=6.0)
            for _ in self.els[e].items
        ]

    def sets(self):
        return 2

    def set(self, k):
        if k == 0:
            return _ns(name="fixed", type="node", items=[1, 4])
        return _ns(name="tri", type="element", items=[101])


def test_marc_t16_route(tmp_path, capsys):
    script = _load("marc_t16/t16_to_vtu.py", "t16_to_vtu_test")
    out = tmp_path / "job.pvd"
    assert script.export(_PyPost(), str(out)) == 2
    assert "9999" in capsys.readouterr().out
    steps = list(meshioplusplus.read_sequence(out))
    assert [t for t, _ in steps] == [0.5, 1.0]
    mesh = steps[1][1]
    assert [c.type for c in mesh.cells] == ["quad", "triangle"]
    np.testing.assert_array_equal(mesh.point_data["Temperature"], [20, 21, 22, 23, 24])
    np.testing.assert_array_equal(
        mesh.point_data["Displacement"][:, 0], [0, 2, 4, 6, 8]
    )
    assert [a[0] for a in mesh.cell_data["Equivalent Von Mises Stress"]] == [
        2.5,
        10 / 3,
    ]
    np.testing.assert_array_equal(mesh.cell_data["Stress"][1][0], [1, 2, 3, 4, 5, 6])
    np.testing.assert_array_equal(mesh.point_data["set:fixed"], [1, 0, 0, 1, 0])
    assert [list(a) for a in mesh.cell_data["set:tri"]] == [[0], [1]]
    assert [list(a) for a in mesh.cell_data["marc:element_id"]] == [[100], [101]]


# --------------------------------------------------------------------------- #
# Ansys results through PyDPF                                                 #
# --------------------------------------------------------------------------- #


def _dpf_stub(pv):
    locations = _ns(
        nodal="Nodal",
        elemental="Elemental",
        elemental_nodal="ElementalNodal",
        time_freq="TimeFreq_steps",
    )
    grid = pv.UnstructuredGrid(
        np.array([4, 0, 1, 2, 3, 3, 1, 4, 2]),
        np.array([pv.CellType.TETRA, pv.CellType.TRIANGLE]),
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0]], float),
    )
    node_ids = [11, 12, 13, 14, 15]
    element_ids = [7, 9]

    class Field:
        def __init__(self, location, ids, data):
            self.location = location
            self.scoping = _ns(ids=ids)
            self.data = np.asarray(data, dtype=float)

    class Container:
        def __init__(self, by_time):
            self.by_time = by_time

        def get_fields(self, label):
            return self.by_time.get(label["time"], [])

    class Operator:
        def __init__(self, name):
            self.name = name
            self.pins = {}

        def connect(self, pin, value):
            self.pins[pin] = value

        def get_output(self, pin, kind):
            sets = self.pins[0].ids
            if self.name == "U":
                return Container(
                    {
                        s: [Field("Nodal", node_ids[::-1], np.full((5, 3), s))]
                        for s in sets
                    }
                )
            # Stress, elemental-nodal: averaged to the nodes when pin 9 asks.
            assert self.pins.get(9) == "Nodal"
            return Container(
                {s: [Field("Nodal", node_ids, np.full((5, 6), 10.0 * s))] for s in sets}
            )

    results = [
        _ns(name="displacement", operator_name="U", native_location="Nodal"),
        _ns(name="stress", operator_name="S", native_location="ElementalNodal"),
    ]

    class Model:
        def __init__(self, path):
            self.metadata = _ns(
                meshed_region=_ns(
                    grid=grid,
                    nodes=_ns(scoping=_ns(ids=node_ids)),
                    elements=_ns(scoping=_ns(ids=element_ids)),
                ),
                time_freq_support=_ns(time_frequencies=_ns(data=[0.25, 1.0])),
                result_info=_ns(available_results=results),
            )

        def operator(self, name):
            return Operator(name)

    return _ns(
        Model=Model,
        locations=locations,
        Scoping=lambda ids, location: _ns(ids=ids, location=location),
        types=_ns(fields_container="fields_container"),
    )


def test_ansys_dpf_route(tmp_path):
    pv = pytest.importorskip("pyvista")
    pytest.importorskip("h5py")
    script = _load("ansys_dpf/rst_to_vtkhdf.py", "rst_to_vtkhdf_test")
    dpf = _dpf_stub(pv)
    out = tmp_path / "file.vtkhdf"
    meshioplusplus.write_sequence(
        out, script.steps(dpf.Model("file.rst"), dpf, ["displacement", "stress"])
    )
    steps = list(meshioplusplus.read_sequence(out))
    assert [t for t, _ in steps] == [0.25, 1.0]
    mesh = steps[1][1]
    np.testing.assert_array_equal(mesh.point_data["displacement"], np.full((5, 3), 2))
    np.testing.assert_array_equal(mesh.point_data["stress"], np.full((5, 6), 20.0))
    np.testing.assert_array_equal(mesh.point_data["dpf:node_id"], [11, 12, 13, 14, 15])


# --------------------------------------------------------------------------- #
# Femap .modfem and Tecplot .szplt                                            #
# --------------------------------------------------------------------------- #


def test_femap_route_calls():
    script = _load("femap/export_neutral.py", "export_neutral_test")
    calls = []

    class Femap:
        def feFileOpen(self, skip_save, path):
            calls.append(("open", skip_save, path))
            return script.FE_OK

        def feFileWriteNeutral(self, *args):
            calls.append(("write",) + args)
            return script.FE_OK

    script.export(Femap(), "model.modfem", "model.neu", "2401")
    assert calls[0][0] == "open" and calls[0][2].endswith("model.modfem")
    write = calls[1]
    assert write[1] == 0 and write[2].endswith("model.neu") and write[8] == 2401.0

    class Failing(Femap):
        def feFileOpen(self, skip_save, path):
            return 0

    with pytest.raises(SystemExit, match="could not open"):
        script.export(Failing(), "model.modfem", "model.neu", "2401")


def test_tecplot_route_calls(tmp_path, monkeypatch):
    script = _load("tecplot/szplt_to_plt.py", "szplt_to_plt_test")
    saved = []
    tp = _ns(
        new_layout=lambda: None,
        data=_ns(
            load_tecplot_szl=lambda path: ("dataset", path),
            save_tecplot_plt=lambda out, dataset: saved.append((out, dataset)),
        ),
    )
    out = script.convert(tp, str(tmp_path / "run.szplt"))
    assert out == str(tmp_path / "run.plt")
    assert saved == [(out, ("dataset", str(tmp_path / "run.szplt")))]
    monkeypatch.setitem(sys.modules, "tecplot", tp)
    script.main([str(tmp_path / "a.szplt")])
    assert saved[-1][0] == str(tmp_path / "a.plt")
