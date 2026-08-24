"""
I/O for the COMSOL ``.mphtxt`` text mesh format, following FEconv
<https://github.com/victorsndvg/FEconv>.

The file stores a version, tag/type name tables, then one or more mesh objects.
A mesh object holds the space dimension, the node coordinates and a series of
element-type blocks (``tet``, ``hex``, ``tri2`` …), each with its connectivity
and a per-element "geometric entity index" (exposed as
``cell_data["mphtxt:geom"]``).  Comments run from ``#`` to end of line.
"""

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._files import open_file
from .._mesh import CellBlock, Mesh
from .._provenance import TAG as _PROVENANCE_TAG

_comsol_to_meshio_type = {
    "vtx": "vertex",
    "edg": "line",
    "tri": "triangle",
    "quad": "quad",
    "tet": "tetra",
    "prism": "wedge",
    "pyr": "pyramid",
    "hex": "hexahedron",
    "edg2": "line3",
    "tri2": "triangle6",
    "quad2": "quad9",
    "tet2": "tetra10",
    "prism2": "wedge18",
    "hex2": "hexahedron27",
}
_meshio_to_comsol_type = {v: k for k, v in _comsol_to_meshio_type.items()}

# COMSOL node ordering -> meshio (self-inverse swaps); identity otherwise.
_perm = {
    "quad": [0, 1, 3, 2],
    "hexahedron": [0, 1, 3, 2, 4, 5, 7, 6],
}


class _Cursor:
    def __init__(self, tokens):
        self.t = tokens
        self.i = 0

    def tok(self):
        v = self.t[self.i]
        self.i += 1
        return v

    def integer(self):
        return int(self.tok())

    def string(self):
        # length-prefixed name; names are single tokens in practice
        self.integer()
        return self.tok()


def read(filename):
    with open_file(filename, "r") as f:
        tokens = []
        for line in f:
            tokens.extend(line.split("#", 1)[0].split())
    c = _Cursor(tokens)

    c.integer()  # version major
    c.integer()  # version minor
    for _ in range(c.integer()):  # tags
        c.string()
    n_types = c.integer()
    for _ in range(n_types):  # type names
        c.string()

    points = np.empty((0, 3))
    cells = []
    cell_geom = []

    for _ in range(n_types):
        c.integer(), c.integer(), c.integer()  # object type indices
        c.string()  # class name ("Mesh")
        c.integer()  # object version
        sdim = c.integer()
        n_points = c.integer()
        lowest = c.integer()
        points = np.array([float(c.tok()) for _ in range(n_points * sdim)]).reshape(
            n_points, sdim
        )

        for _ in range(c.integer()):  # element types
            comsol_type = c.string()
            if comsol_type not in _comsol_to_meshio_type:
                raise ReadError(f"mphtxt: unknown element type '{comsol_type}'")
            mtype = _comsol_to_meshio_type[comsol_type]
            n_nodes = c.integer()
            n_elem = c.integer()
            conn = np.array(
                [c.integer() for _ in range(n_elem * n_nodes)], dtype=int
            ).reshape(n_elem, n_nodes)
            conn = conn - lowest
            p = _perm.get(mtype)
            if p is not None:
                conn = conn[:, p]
            n_par_per = c.integer()
            n_par = c.integer()
            for _ in range(n_par * n_par_per):
                c.tok()
            n_geom = c.integer()
            gvals = np.array([c.integer() for _ in range(n_geom)], dtype=int)
            n_ud = c.integer()
            for _ in range(n_ud * 2):
                c.integer()
            cells.append(CellBlock(mtype, conn))
            cell_geom.append(gvals)
        break  # only the first mesh object

    cell_data = {"mphtxt:geom": cell_geom} if cell_geom else {}
    return Mesh(points, cells, cell_data=cell_data)


def write(filename, mesh):
    sdim = mesh.points.shape[1]
    blocks = []
    for k, cb in enumerate(mesh.cells):
        if cb.type not in _meshio_to_comsol_type:
            warn(f"mphtxt does not support '{cb.type}' cells. Skipping.")
            continue
        blocks.append((k, cb))

    geom_data = mesh.cell_data.get("mphtxt:geom")

    with open_file(filename, "w") as f:
        f.write(f"# {_PROVENANCE_TAG}\n\n")
        f.write("0 1\n")  # version
        f.write("1 # number of tags\n5 mesh1\n")
        f.write("1 # number of types\n3 obj\n\n")

        # object
        f.write("0 0 1\n")
        f.write("4 Mesh # class\n")
        f.write("2 # version\n")
        f.write(f"{sdim} # sdim\n")
        f.write(f"{len(mesh.points)} # number of mesh points\n")
        f.write("1 # lowest mesh point index\n\n")
        f.write("# Mesh point coordinates\n")
        for pt in mesh.points:
            f.write(" ".join(repr(float(x)) for x in pt) + "\n")
        f.write("\n")

        f.write(f"{len(blocks)} # number of element types\n\n")
        for ti, (k, cb) in enumerate(blocks):
            comsol_type = _meshio_to_comsol_type[cb.type]
            data = cb.data
            p = _perm.get(cb.type)
            if p is not None:
                data = data[:, p]
            n_nodes = data.shape[1]
            n_elem = data.shape[0]
            f.write(f"# Type #{ti + 1}\n\n")
            f.write(f"{len(comsol_type)} {comsol_type} # type name\n\n")
            f.write(f"{n_nodes} # number of nodes per element\n")
            f.write(f"{n_elem} # number of elements\n")
            f.write("# Elements\n")
            for row in data + 1:  # lowest index is 1
                f.write(" ".join(str(v) for v in row) + "\n")
            f.write(f"\n{n_nodes} # number of parameter values per element\n")
            f.write("0 # number of parameters\n# Parameters\n\n")
            if geom_data is not None and k < len(geom_data):
                gvals = np.asarray(geom_data[k], dtype=int)
            else:
                gvals = np.zeros(n_elem, dtype=int)
            f.write(f"{len(gvals)} # number of geometric entity indices\n")
            f.write("# Geometric entity indices\n")
            for g in gvals:
                f.write(f"{int(g)}\n")
            f.write("\n0 # number of up/down pairs\n# Up/down\n\n")
