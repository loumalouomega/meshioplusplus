"""Which format a write picks when none is named.

A write resolves its format from the path alone, among the formats that can
write: a read-only format registered for the same extension must never shadow a
writer (until v16.17.0 ``write("a.dat", mesh)`` failed with "Unknown format
'marc'", Marc's read-only ``.dat`` registration standing in front of Tecplot's),
and an existing file's content never chooses the format of the file that
replaces it (natively, ``resolve_write_format``).
"""

import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._helpers import (
    _write_format_for_path,
    _writer_map,
    extension_to_filetypes,
)

# Native registry keys that Python spells differently.
_NATIVE_ALIASES = {"ansysinp": "ansysInp", "dolfin": "dolfin-xml"}
# Python picks among ".msh"'s writers by the mesh's tags (_pick_best_format);
# with no mesh in hand it takes the first, and the native default is gmsh.
_KNOWN_DIVERGENCE = {".msh"}


def _tri():
    return meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2]]))],
    )


@pytest.mark.parametrize("ext", sorted(extension_to_filetypes))
def test_a_write_picks_a_writer_when_the_extension_has_one(ext):
    candidates = extension_to_filetypes[ext]
    chosen = _write_format_for_path(pathlib.Path("x" + ext))
    if any(c in _writer_map for c in candidates):
        assert chosen in _writer_map, (ext, candidates)


def test_python_and_native_extension_defaults_agree():
    native = _core.registry_formats()["extensions"]
    for ext, candidates in sorted(extension_to_filetypes.items()):
        if ext in _KNOWN_DIVERGENCE or ext not in native:
            continue
        if not any(c in _writer_map for c in candidates):
            continue
        want = _NATIVE_ALIASES.get(native[ext], native[ext])
        assert _write_format_for_path(pathlib.Path("x" + ext)) == want, ext


def test_dat_writes_tecplot(tmp_path):
    out = tmp_path / "a.dat"
    meshioplusplus.write(out, _tri())
    assert meshioplusplus.read(out, file_format="tecplot").points.shape == (3, 3)
    assert _core.resolve_write_format(str(out), "") == "tecplot"


def test_an_existing_file_never_chooses_the_written_format(tmp_path):
    # A Marc deck at the path still resolves to Marc for reading, but a write
    # to the same path is Tecplot's, on both surfaces.
    deck = tmp_path / "deck.dat"
    deck.write_text(
        "title\nextended\nend\nconnectivity\n"
        + "         1         0         1\n"
        + "         1         6         1         2         3\n"
        + "coordinates\n         3         3         0         1\n"
        + "         1"
        + "".join(f"{v:20.12e}" for v in (0.0, 0.0, 0.0))
        + "\n"
        + "         2"
        + "".join(f"{v:20.12e}" for v in (1.0, 0.0, 0.0))
        + "\n"
        + "         3"
        + "".join(f"{v:20.12e}" for v in (0.0, 1.0, 0.0))
        + "\n"
        + "end option\n"
    )
    assert _core.registry_formats()["extensions"][".dat"] == "tecplot"
    assert meshioplusplus.read(deck).cells[0].type == "triangle"  # read as Marc
    assert _core.resolve_write_format(str(deck), "") == "tecplot"
    meshioplusplus.write(deck, _tri())
    assert meshioplusplus.read(deck, file_format="tecplot").points.shape == (3, 3)


def test_naming_a_read_only_format_says_so(tmp_path):
    with pytest.raises(meshioplusplus.WriteError, match="can be read but not written"):
        meshioplusplus.write(tmp_path / "a.frd", _tri())
