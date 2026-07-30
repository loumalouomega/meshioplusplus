"""
I/O for DOLFIN's XML format, cf.
<https://people.sc.fsu.edu/~jburkardt/data/dolfin_xml/dolfin_xml.html>.
"""

import os
import pathlib
import re
from xml.etree import ElementTree as ET

import numpy as np

from .._common import warn
from .._exceptions import ReadError, WriteError
from .._mesh import Mesh


def _read_mesh(filename):
    dolfin_to_meshio_type = {"triangle": ("triangle", 3), "tetrahedron": ("tetra", 4)}

    # Use iterparse() to avoid loading the entire file via parse(). iterparse()
    # allows to discard elements (via clear()) after they have been processed.
    # See <https://stackoverflow.com/a/326541/353337>.
    dim = None
    points = None
    keys = None
    cell_type = None
    num_nodes_per_cell = None
    cells = None
    cell_tags = None
    for event, elem in ET.iterparse(filename, events=("start", "end")):
        if event == "end":
            continue

        if elem.tag == "dolfin":
            # Don't be too strict with the assertion. Some mesh files don't have the
            # proper tags.
            # assert elem.attrib['nsmap'] \
            #     == '{\'dolfin\': \'https://fenicsproject.org/\'}'
            pass
        elif elem.tag == "mesh":
            dim = int(elem.attrib["dim"])
            cell_type, num_nodes_per_cell = dolfin_to_meshio_type[
                elem.attrib["celltype"]
            ]
            cell_tags = [f"v{i}" for i in range(num_nodes_per_cell)]
        elif elem.tag == "vertices":
            if dim is None:
                raise ReadError("Expected `mesh` before `vertices`")
            points = np.empty((int(elem.attrib["size"]), dim))
            keys = ["x", "y"]
            if dim == 3:
                keys += ["z"]
        elif elem.tag == "vertex":
            if points is None or keys is None:
                raise ReadError("Expected `vertices` before `vertex`")
            k = int(elem.attrib["index"])
            points[k] = [elem.attrib[key] for key in keys]
        elif elem.tag == "cells":
            if cell_type is None or num_nodes_per_cell is None:
                raise ReadError("Expected `mesh` before `cells`")
            cells = [
                (
                    cell_type,
                    np.empty((int(elem.attrib["size"]), num_nodes_per_cell), dtype=int),
                )
            ]
        elif elem.tag in ["triangle", "tetrahedron"]:
            k = int(elem.attrib["index"])
            assert cells is not None
            assert cell_tags is not None
            cells[0][1][k] = [elem.attrib[t] for t in cell_tags]
        else:
            warn(f"Unknown entry {elem.tag}. Ignoring.")

        elem.clear()

    return points, cells, cell_type


def _read_mesh_functions(filename):
    """Read the sibling ``<stem>_<name>.xml`` mesh functions.

    Returns ``(point_data, cell_data)``. A ``mesh_function``'s ``dim`` attribute
    is the topological dimension of the entities it is defined on, so ``dim="0"``
    means *vertices* and anything else means cells — which is the whole
    discriminator, and why point data needs no new file convention. Twin of the
    corresponding loop in ``src/cpp/src/formats/dolfin.cpp``.
    """
    dolfin_type_to_numpy_type = {
        "int": np.dtype("int"),
        "float": np.dtype("float"),
        "uint": np.dtype("uint"),
    }

    point_data = {}
    cell_data = {}
    dir_name = pathlib.Path(filename).resolve().parent

    # Loop over all files in the same directory as `filename`.
    basename = pathlib.Path(filename).stem
    for f in os.listdir(dir_name):
        # Check if there are files by the name "<filename>_*.xml"; if yes,
        # extract the * pattern and make it the name of the data set.
        out = re.match(f"{basename}_([^\\.]+)\\.xml", f)
        if not out:
            continue
        name = out.group(1)

        parser = ET.XMLParser()
        tree = ET.parse((dir_name / f).as_posix(), parser)
        root = tree.getroot()

        mesh_functions = list(root)
        if len(mesh_functions) != 1:
            raise ReadError("Can only handle one mesh function")
        mesh_function = mesh_functions[0]

        if mesh_function.tag != "mesh_function":
            raise ReadError()
        size = int(mesh_function.attrib["size"])
        dtype = dolfin_type_to_numpy_type[mesh_function.attrib["type"]]
        data = np.empty(size, dtype=dtype)
        for child in mesh_function:
            if child.tag != "entity":
                raise ReadError()
            idx = int(child.attrib["index"])
            data[idx] = child.attrib["value"]

        if int(mesh_function.attrib.get("dim", -1)) == 0:
            point_data[name] = data
        else:
            if name not in cell_data:
                cell_data[name] = []
            cell_data[name].append(data)

    return point_data, cell_data


def read(filename):
    points, cells, _ = _read_mesh(filename)
    point_data, cell_data = _read_mesh_functions(filename)
    return Mesh(points, cells, point_data=point_data, cell_data=cell_data)


def _write_mesh(filename, points, cell_type, cells):
    stripped_cells = [c for c in cells if c.type == cell_type]

    meshio_to_dolfin_type = {"triangle": "triangle", "tetra": "tetrahedron"}

    if any(c.type != cell_type for c in cells):
        discarded_cell_types = {c.type for c in cells if c.type != cell_type}
        warn(
            "DOLFIN XML can only handle one cell type at a time. "
            + f"Using {cell_type}, discarding {', '.join(discarded_cell_types)}.",
        )

    dim = points.shape[1]
    if dim not in [2, 3]:
        raise WriteError(f"Can only write dimension 2, 3, got {dim}.")

    coord_names = ["x", "y"]
    if dim == 3:
        coord_names += ["z"]

    with open(filename, "w") as f:
        f.write("<dolfin nsmap=\"{'dolfin': 'https://fenicsproject.org/'}\">\n")
        ct = meshio_to_dolfin_type[cell_type]
        f.write(f'  <mesh celltype="{ct}" dim="{dim}">\n')

        num_points = len(points)
        f.write(f'    <vertices size="{num_points}">\n')
        for idx, point in enumerate(points):
            s = " ".join(f'{xyz}="{p}"' for xyz, p in zip("xyz", point))
            f.write(f'      <vertex index="{idx}" {s} />\n')
        f.write("    </vertices>\n")

        num_cells = 0
        for c in stripped_cells:
            num_cells += len(c.data)

        f.write(f'    <cells size="{num_cells}">\n')
        idx = 0
        for cell_block in stripped_cells:
            type_string = meshio_to_dolfin_type[cell_block.type]
            for cell in cell_block.data:
                s = " ".join(f'v{k}="{c}"' for k, c in enumerate(cell))
                f.write(f'      <{type_string} index="{idx}" {s} />\n')
                idx += 1
        f.write("    </cells>\n")
        f.write("  </mesh>\n")
        f.write("</dolfin>")


def _numpy_type_to_dolfin_type(dtype):
    types = {
        "int": [np.int8, np.int16, np.int32, np.int64],
        "uint": [np.uint8, np.uint16, np.uint32, np.uint64],
        "float": [np.float16, np.float32, np.float64],
    }
    for key, numpy_types in types.items():
        for numpy_type in numpy_types:
            if np.issubdtype(dtype, numpy_type):
                return key

    raise WriteError("Could not convert NumPy data type to DOLFIN data type.")


def _write_mesh_function(filename, dim, values):
    dolfin = ET.Element("dolfin", nsmap={"dolfin": "https://fenicsproject.org/"})

    mesh_function = ET.SubElement(
        dolfin,
        "mesh_function",
        type=_numpy_type_to_dolfin_type(values.dtype),
        dim=str(dim),
        size=str(len(values)),
    )

    for k, value in enumerate(values):
        ET.SubElement(mesh_function, "entity", index=str(k), value=str(value))

    tree = ET.ElementTree(dolfin)
    tree.write(filename)


def write(filename, mesh):
    warn("DOLFIN XML is a legacy format. Consider using XDMF instead.")

    if any("tetra" == c.type for c in mesh.cells):
        cell_type = "tetra"
    elif any("triangle" == c.type for c in mesh.cells):
        cell_type = "triangle"
    else:
        raise WriteError(
            "DOLFIN XML only supports triangles and tetrahedra. "
            "Consider using XDMF instead."
        )

    _write_mesh(filename, mesh.points, cell_type, mesh.cells)

    fname = os.path.splitext(filename)[0]
    dim = 2 if mesh.points.shape[1] == 2 or all(mesh.points[:, 2] == 0) else 3

    for name, lst in mesh.cell_data.items():
        for data in lst:
            _write_mesh_function(f"{fname}_{name}.xml", dim, np.array(data))

    # Point data, as `dim="0"` mesh functions -- vertices are the topological
    # entities of dimension 0, so this is the format's own notion rather than a
    # meshio++ convention. A name used by *both* locations would want the same
    # sibling file, and cell data has always owned it, so the point array is
    # skipped with a warning rather than silently clobbering it.
    for name, data in mesh.point_data.items():
        if name in mesh.cell_data:
            warn(
                f"DOLFIN: point_data '{name}' collides with a cell_data array of "
                "the same name (both want the same sibling file); not written."
            )
            continue
        values = np.asarray(data)
        if values.ndim != 1:
            warn(
                f"DOLFIN: point_data '{name}' has {values.shape[1:]} components; "
                "a mesh function is scalar per entity, so it is not written."
            )
            continue
        _write_mesh_function(f"{fname}_{name}.xml", 0, values)
