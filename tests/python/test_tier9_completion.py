"""Roadmap §9's last pieces, across both engines where a format has two: the
CRS provider (PROJ) against pyproj, ``meshio:unit:<array>`` through Selafin,
MED and CGNS (SIDS ``DimensionalExponents``), ``meshio:crs`` through LAS and
UGRID netCDF, the ``SurfaceQuery`` handle (9.3.2), the hexahedron mixed-green
closure (9.6), and the compiled-in graph partitioners and boolean provider."""

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _units

_core = pytest.importorskip("meshioplusplus._core")

# ---- units ---------------------------------------------------------------


def _tri_mesh():
    pts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    mesh = mio.Mesh(pts, [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))])
    mesh.point_data["p"] = np.array([1.0, 2.0, 3.0, 4.0])
    mio.set_field_unit(mesh, "p", "MPa")
    return mesh


def _engines(module, core_read, core_write):
    py = __import__(f"meshioplusplus.{module}", fromlist=["read", "write"])
    return {
        "core": (getattr(_core, core_read), getattr(_core, core_write)),
        "python": (py.read, py.write),
    }


@pytest.mark.parametrize("writer", ["core", "python"])
@pytest.mark.parametrize("reader", ["core", "python"])
def test_selafin_carries_the_field_unit(tmp_path, writer, reader):
    eng = _engines("selafin._selafin", "selafin_read", "selafin_write")
    path = str(tmp_path / "u.slf")
    eng[writer][1](path, _tri_mesh())
    back = eng[reader][0](path)
    unit = mio.field_unit(back, "p")
    assert unit is not None and unit[0] == pytest.approx(1e6)
    assert unit[1] == (-1.0, 1.0, -2.0, 0.0, 0.0, 0.0, 0.0)
    raw = np.asarray(back.field_data["selafin:unit:p"], np.uint8).tobytes()
    assert raw.decode("latin-1").strip() == "MPa"


def test_selafin_reads_telemac_spellings():
    np.testing.assert_array_equal(
        _units.field_unit_from_text("M/S"), [1.0, 1, 0, -1, 0, 0, 0, 0]
    )
    np.testing.assert_array_equal(
        _units.field_unit_from_text("KG/M3"), [1.0, -3, 1, 0, 0, 0, 0, 0]
    )


@pytest.mark.parametrize("writer", ["core", "python"])
def test_med_writes_the_field_unit_and_reads_it_back(tmp_path, writer):
    pytest.importorskip("h5py")
    from meshioplusplus.med import _med

    mesh = _tri_mesh()
    path = str(tmp_path / "u.med")
    if writer == "core":
        mio.write(path, mesh)
    else:
        _med.write(path, mesh)
    back = _med.read(path)
    unit = mio.field_unit(back, "p")
    assert unit is not None and unit[0] == pytest.approx(1e6)


@pytest.mark.parametrize("writer", ["core", "python"])
@pytest.mark.parametrize("reader", ["core", "python"])
def test_cgns_states_the_unit_as_sids_exponents(tmp_path, writer, reader):
    h5py = pytest.importorskip("h5py")
    from meshioplusplus.cgns import _cgns

    mesh = _tri_mesh()
    mio.set_field_unit(mesh, "p", "mm/s")
    path = str(tmp_path / "u.cgns")
    if writer == "core":
        _core.cgns_write(path, mesh, 4)
    else:
        _cgns.write(path, mesh)
    with h5py.File(path, "r") as f:
        node = f["Base/Zone1/FlowSolution/p"]
        np.testing.assert_array_equal(
            node["DimensionalExponents/ data"][()], [0.0, 1.0, -1.0, 0.0, 0.0]
        )
        np.testing.assert_allclose(node["DataConversion/ data"][()], [1e-3, 0.0])
        assert "meshio:unit:p" not in f["Base/Zone1"].get("FieldData", {})
    back = _core.cgns_read(path) if reader == "core" else _cgns.read(path)
    unit = mio.field_unit(back, "p")
    assert unit is not None and unit[0] == pytest.approx(1e-3)
    assert unit[1] == (1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0)


# ---- CRS -----------------------------------------------------------------

_WKT = (
    'GEOGCRS["WGS 84",DATUM["World Geodetic System 1984",ELLIPSOID["WGS 84",'
    '6378137,298.257223563]],CS[ellipsoidal,2],AXIS["latitude",north],'
    'AXIS["longitude",east],UNIT["degree",0.0174532925199433],ID["EPSG",4326]]'
)


def _cloud(crs):
    pts = np.array([[2.35, 48.85, 30.0], [2.36, 48.86, 31.0]])
    mesh = mio.Mesh(pts, [("vertex", np.array([[0], [1]]))])
    mio.set_mesh_crs(mesh, crs)
    return mesh


@pytest.mark.parametrize("writer", ["core", "python"])
@pytest.mark.parametrize("reader", ["core", "python"])
def test_las_carries_a_wkt_crs(tmp_path, writer, reader):
    eng = _engines("las._las", "las_read", "las_write")
    path = str(tmp_path / "c.las")
    eng[writer][1](path, _cloud(_WKT))
    assert mio.mesh_crs(eng[reader][0](path)) == _WKT


def test_las_warns_on_an_epsg_code_it_cannot_write(tmp_path, capsys):
    from meshioplusplus.las import _las

    path = str(tmp_path / "c.las")
    _las.write(path, _cloud("EPSG:4326"))
    assert "not WKT" in capsys.readouterr().err.replace("\n", " ")
    assert mio.mesh_crs(_las.read(path)) is None


@pytest.mark.parametrize("crs", ["EPSG:3857", _WKT])
@pytest.mark.parametrize("writer", ["core", "python"])
@pytest.mark.parametrize("reader", ["core", "python"])
def test_ugrid_nc_carries_the_crs(tmp_path, crs, writer, reader):
    pytest.importorskip("netCDF4")
    if not hasattr(_core, "ugrid_nc_write"):
        pytest.skip("needs the netCDF-enabled core")
    eng = _engines("ugrid_nc._ugrid_nc", "ugrid_nc_read", "ugrid_nc_write")
    mesh = _tri_mesh()
    mesh.point_data.clear()
    mesh.field_data.clear()
    mio.set_mesh_crs(mesh, crs)
    path = str(tmp_path / "c.nc")
    eng[writer][1](path, mesh)
    assert mio.mesh_crs(eng[reader][0](path)) == crs


@pytest.mark.skipif(not mio.has_crs_provider(), reason="needs the PROJ provider")
def test_the_proj_provider_matches_pyproj():
    pyproj = pytest.importorskip("pyproj")
    mesh = _cloud("EPSG:4326")
    out = mio.convert_coordinates(mesh, "crs", target_crs="EPSG:3857")
    t = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    x, y, z = t.transform(mesh.points[:, 0], mesh.points[:, 1], mesh.points[:, 2])
    np.testing.assert_allclose(
        out.points, np.column_stack([x, y, z]), rtol=0, atol=1e-6
    )
    assert mio.mesh_crs(out) == "EPSG:3857"


# ---- SurfaceQuery (9.3.2) ------------------------------------------------


def test_surface_query_matches_cast_rays_bit_for_bit():
    sphere = mio.icosa_sphere()
    rng = np.random.default_rng(3)
    origins = rng.normal(size=(50, 3)) * 0.1
    dirs = rng.normal(size=(50, 3))
    query = mio.SurfaceQuery(sphere)
    assert query.num_triangles == len(sphere.cells[0].data)
    a = mio.cast_rays(sphere, origins, dirs)
    b = query.cast_rays(origins, dirs)
    for key in ("hit", "distance", "point", "normal", "cell"):
        np.testing.assert_array_equal(a[key], b[key])
    assert a["num_hits"] == b["num_hits"] == 50
    c = mio.cast_rays(sphere, origins, dirs, all_hits=True)
    d = query.cast_rays(origins, dirs, all_hits=True)
    np.testing.assert_array_equal(c["ray_offset"], d["ray_offset"])


def test_surface_query_closest_points():
    sphere = mio.icosa_sphere()
    query = mio.SurfaceQuery(sphere)
    pts = np.array([[0.0, 0.0, 3.0], [0.0, 0.0, 0.0]])
    res = query.closest_points(pts)
    assert res["distance"][0] == pytest.approx(2.0, abs=0.1)
    assert res["distance"][1] == pytest.approx(1.0, abs=0.1)
    assert np.all(res["cell"] >= 0)


def test_surface_query_refuses_a_mesh_without_a_surface():
    with pytest.raises(ValueError):
        mio.SurfaceQuery(mio.box(dims=(1, 1, 1)))


# ---- mixed green on hexahedra (9.6) --------------------------------------


@pytest.mark.parametrize("cells", [[13], [0], [13, 14], [0, 4, 13]])
def test_hexahedron_mixed_green_is_conforming_and_keeps_the_volume(cells):
    grid = mio.box(dims=(3, 3, 3))
    grid.regions = []
    out = mio.refine(grid, cells=cells, closure="mixed_green")
    types = {cb.type for cb in out.cells}
    assert "hexahedron" in types and "pyramid" in types
    _, report = mio.check_conformity(out, return_report=True)
    assert report["conforming"] and report["num_hanging_nodes"] == 0
    assert mio.mass_properties(out)["total"]["measure"] == pytest.approx(1.0)
    pts = out.points
    for cb in out.cells:
        c = cb.data
        if cb.type == "tetra":
            vol = np.einsum(
                "ij,ij->i",
                np.cross(pts[c[:, 1]] - pts[c[:, 0]], pts[c[:, 2]] - pts[c[:, 0]]),
                pts[c[:, 3]] - pts[c[:, 0]],
            )
            assert vol.min() > 0


# ---- providers -----------------------------------------------------------


@pytest.mark.parametrize("name", ["metis", "scotch"])
def test_a_compiled_in_graph_partitioner_partitions(name):
    if name not in mio.graph_partitioner_names():
        pytest.skip(f"the core was built without {name}")
    grid = mio.box(dims=(4, 4, 4))
    grid.regions = []
    labels = mio.partition_labels(grid, 4, method=name)
    flat = np.concatenate([np.asarray(b).ravel() for b in labels])
    assert sorted(set(flat.tolist())) == [0, 1, 2, 3]


@pytest.mark.skipif(not mio.has_boolean_provider(), reason="needs a boolean provider")
@pytest.mark.parametrize("op", ["union", "difference", "intersection"])
def test_mesh_boolean_volumes_add_up(op):
    a = mio.icosa_sphere()
    b = mio.Mesh(a.points + [0.5, 0.0, 0.0], a.cells)
    va = mio.mass_properties(a)["total"]["measure"]
    vol = {
        o: mio.mass_properties(mio.mesh_boolean(a, b, o))["total"]["measure"]
        for o in ("union", "difference", "intersection")
    }
    # |A u B| = |A| + |B| - |A n B| and |A - B| = |A| - |A n B|.
    assert vol["union"] == pytest.approx(2 * va - vol["intersection"], rel=1e-9)
    assert vol["difference"] == pytest.approx(va - vol["intersection"], rel=1e-9)
    assert vol[op] > 0
