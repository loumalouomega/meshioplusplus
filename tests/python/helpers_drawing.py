"""Shared pieces of the SVG/TikZ reader tests.

Both readers have a C++ core and a Python twin that must agree exactly, so
every test reads a file through both and compares them (:func:`read_both`)
before looking at the mesh itself.
"""

import numpy as np

import meshioplusplus

# A marked SVG document around a body (the marker a planar write leaves).
SVG_MARKER = (
    "<!--\nWritten by meshio++ v1.0.0\nDrawing: planar scale=1 offset=0,10\n-->"
)


def marked_svg(body, root_extra=""):
    return (
        '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" '
        f'viewBox="0 0 10 10"{root_extra}>{SVG_MARKER}{body}</svg>'
    )


def marked_tikz(body):
    return (
        "% Written by meshio++ v1.0.0\n% Drawing: planar\n"
        f"\\begin{{tikzpicture}}\n{body}\n\\end{{tikzpicture}}\n"
    )


def readers(fmt):
    """``(core, python)`` read callables of ``fmt`` ("svg" or "tikz")."""
    twin = getattr(meshioplusplus, fmt)
    module = getattr(twin, "_" + fmt)
    return getattr(meshioplusplus._core, fmt + "_read"), module.read


def assert_same(a, b):
    """The two meshes are identical: points, blocks and cell data."""
    assert a.points.dtype == b.points.dtype == np.float64
    assert a.points.shape == b.points.shape
    assert np.array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for ca, cb in zip(a.cells, b.cells):
        assert np.array_equal(np.asarray(ca.data), np.asarray(cb.data))
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for xa, xb in zip(a.cell_data[name], b.cell_data[name]):
            assert xa.dtype == xb.dtype
            assert np.array_equal(xa, xb)
    assert not a.point_data and not b.point_data


def read_both(fmt, path, **kwargs):
    """Reads ``path`` through the C++ core and the Python twin, asserts they agree
    and returns the mesh."""
    core, py = readers(fmt)
    a = core(str(path), **kwargs)
    b = py(str(path), **kwargs)
    assert_same(a, b)
    return a


def cell_coordinates(mesh):
    """Each cell as its sorted tuple of rounded corner coordinates, per type --
    point numbering and corner order do not matter to a drawing."""
    out = {}
    for block in mesh.cells:
        for row in np.asarray(block.data):
            key = tuple(sorted(tuple(np.round(mesh.points[i][:2], 4)) for i in row))
            out.setdefault(block.type, []).append(key)
    return {t: sorted(v) for t, v in out.items()}
