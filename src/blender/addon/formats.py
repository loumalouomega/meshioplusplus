"""Format lists for the file dialog, derived from meshio++'s own registry.

Everything here comes from the public :func:`meshioplusplus.formats`. Nothing
reaches into ``_helpers.reader_map``: that is what
``tools/paraview-meshioplusplus-plugin.py`` does, and a rename there breaks it
silently at load time inside the host application, where nobody is watching.
"""

from __future__ import annotations

import meshioplusplus as mio

#: Blender stores ``FileSelectParams.filter_glob`` in a fixed-size C buffer.
#: The exact size is 255 or 256 depending on the release; the budget below is
#: deliberately one under the smaller value, so the behaviour is the same
#: either way and no extension is ever cut in half.
GLOB_LIMIT = 254

#: Extensions worth guaranteeing a place in the dialog when the full list does
#: not fit. Ordered by how often they actually turn up in FEA work, not
#: alphabetically -- the tail of a truncated glob is invisible, so this is the
#: only part of the list a user can rely on seeing.
PRIORITY = (
    ".msh",
    ".vtu",
    ".vtk",
    ".stl",
    ".obj",
    ".ply",
    ".inp",
    ".med",
    ".cgns",
    ".xdmf",
    ".e",
    ".exo",
    ".mesh",
    ".off",
    ".unv",
    ".mdpa",
    ".vtp",
    ".su2",
    ".bdf",
    ".case",
)


def _fit(extensions, limit):
    """A ``*.a;*.b`` glob within ``limit`` characters, priority extensions first."""
    ordered = [e for e in PRIORITY if e in extensions]
    ordered += [e for e in extensions if e not in ordered]
    out: list = []
    used = 0
    for ext in ordered:
        piece = len(ext) + 1 + (1 if out else 0)  # "*" + ext, plus a ";"
        if used + piece > limit:
            continue
        out.append("*" + ext)
        used += piece
    return ";".join(out)


def import_glob():
    """The file-dialog filter, trimmed to what Blender's buffer can hold.

    Truncation is silent in Blender, so this trims deliberately and the import
    operator offers a ``show_all_files`` toggle that clears the filter — no
    format is ever unreachable, only absent from the default view.
    """
    extensions = sorted(mio.formats()["extensions"])
    full = ";".join("*" + e for e in extensions)
    return full if len(full) <= GLOB_LIMIT else _fit(extensions, GLOB_LIMIT)


def _items(names, auto_label):
    items = [("AUTO", auto_label, "Infer the format from the file name")]
    items += [(name, name, f"Force the '{name}' format") for name in names]
    return items


def read_items():
    return _items(mio.formats()["readable"], "Automatic")


def write_items():
    return _items(mio.formats()["writable"], "From file extension")
