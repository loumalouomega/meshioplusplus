"""``detect``: the ranked format report (roadmap 2.4.9), both engines."""

import numpy as np
import pytest

import meshioplusplus as mio

_TRI = mio.Mesh(
    np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], dtype=float),
    [("triangle", np.array([[0, 1, 2], [1, 3, 2]]))],
)


def _core():
    try:
        from meshioplusplus import _core
    except ImportError:
        return None
    return _core


def test_a_signature_is_certain(tmp_path):
    path = tmp_path / "a.vtu"
    mio.write(path, _TRI)
    found = mio.detect(path)
    assert found[0][:2] == ("vtu", "certain")


def test_gmsh_22_and_41_are_certain_and_rank_the_others_after(tmp_path):
    for version in ("gmsh22", "gmsh"):
        path = tmp_path / f"{version}.msh"
        mio.write(path, _TRI, file_format=version)
        found = mio.detect(path)
        assert found[0][:2] == ("gmsh", "certain"), found
        assert {f for f, _, _ in found} >= {"ansys", "freefem"}
        assert all(c != "certain" for _, c, _ in found[1:])


def test_abaqus_and_calculix_inputs_are_abaqus(tmp_path):
    abaqus = tmp_path / "a.inp"
    mio.write(abaqus, _TRI, file_format="abaqus")
    calculix = tmp_path / "c.inp"
    calculix.write_text(
        "*NODE, NSET=NALL\n1, 0, 0, 0\n2, 1, 0, 0\n3, 0, 1, 0\n"
        "*ELEMENT, TYPE=S3, ELSET=EALL\n1, 1, 2, 3\n"
    )
    for path in (abaqus, calculix):
        found = mio.detect(path)
        assert found[0][0] == "abaqus", found
        assert found[0][1] in ("certain", "likely")


def test_medit_and_mfem_meshes_resolve(tmp_path):
    medit = tmp_path / "a.mesh"
    mio.write(medit, _TRI, file_format="medit")
    assert mio.detect(medit)[0][0] == "medit"
    mfem = tmp_path / "b.mesh"
    mio.write(mfem, _TRI, file_format="mfem")
    assert mio.detect(mfem)[0][0] == "mfem"


def test_vtk_xml_with_an_xml_extension(tmp_path):
    vtu = tmp_path / "a.vtu"
    mio.write(vtu, _TRI)
    path = tmp_path / "a.xml"
    path.write_bytes(vtu.read_bytes())
    assert mio.detect(path)[0][:2] == ("vtu", "certain")


def test_an_unknown_file_is_an_empty_report(tmp_path):
    path = tmp_path / "x.unknownext"
    path.write_text("nothing here\n")
    assert mio.detect(path) == []


def test_a_wrong_extension_is_a_guess(tmp_path):
    path = tmp_path / "a.stl"
    path.write_text("not a mesh at all\n")
    found = mio.detect(path)
    assert found and all(c == "guess" for _, c, _ in found)


@pytest.mark.parametrize("ext", [".msh", ".inp", ".mesh", ".vtu", ".stl", ".dat"])
def test_the_engines_agree_on_core_formats(tmp_path, ext):
    core = _core()
    if core is None:
        pytest.skip("no native core")
    fmt = {".msh": "gmsh", ".inp": "abaqus", ".mesh": "medit", ".dat": "tecplot"}.get(
        ext
    )
    path = tmp_path / f"a{ext}"
    mio.write(path, _TRI, file_format=fmt)
    native = [tuple(e) for e in core.detect_format(str(path))]
    full = mio.detect(path)
    # Python-only formats are candidates of the twin alone.
    py = [e for e in full if e[0] in {n[0] for n in native}]
    assert native == py
    assert native[0] == full[0]


def test_the_cli_reports_and_exits(tmp_path, capsys):
    path = tmp_path / "a.vtu"
    mio.write(path, _TRI)
    assert mio._cli.main(["detect", str(path)]) == 0
    assert "vtu  certain" in capsys.readouterr().out
    other = tmp_path / "x.unknownext"
    other.write_text("?")
    assert mio._cli.main(["detect", str(other)]) == 1
