"""
I/O for XDMF.
https://xdmf.org/index.php/XDMF_Model_and_Format
"""

import os
import pathlib
from io import BytesIO
from xml.etree import ElementTree as ET

import numpy as np

from .._common import cell_data_from_raw, raw_from_cell_data, warn, write_xml
from .._exceptions import ReadError, WriteError
from .._mesh import CellBlock, Mesh, topological_dimension
from .._regions import Region, block_bases
from .common import (
    attribute_type,
    dtype_to_format_string,
    meshio_to_xdmf_type,
    meshio_type_to_xdmf_index,
    numpy_to_xdmf_dtype,
    translate_mixed_cells,
    xdmf_to_meshio_type,
    xdmf_to_numpy_type,
)


def read(filename):
    return XdmfReader(filename).read()


def _block_dim(block):
    """A block's topological dimension (a polyhedron is 3, a polygon 2)."""
    if block.type.startswith("polyhedron"):
        return 3
    if block.type.startswith("polygon"):
        return 2
    return topological_dimension.get(block.type, 3)


_SET_KIND = {"Node": "point", "Cell": "cell", "Face": "side", "Edge": "side"}


def _read_set(reader, c, regions):
    """One ``<Set>`` into ``regions`` (name -> Region), the twin of the core's
    ``xdmf_read_set``: Node/Cell/Face/Edge -> point/cell/side, the cell ids
    first and the cell-local face/edge ids second."""
    name = c.attrib.get("Name", "")
    set_type = c.attrib.get("SetType", "")
    kind = _SET_KIND.get(set_type)
    if kind is None:
        warn(f"xdmf: skipping set '{name}' of SetType '{set_type}'")
        return
    dim, tag = -1, -1
    for info in c.findall("Information"):
        if info.attrib.get("Name") == "meshio++:dim":
            dim = int(info.attrib["Value"])
        elif info.attrib.get("Name") == "meshio++:tag":
            tag = int(info.attrib["Value"])
    if c.find("Attribute") is not None:
        warn(f"xdmf: set '{name}' carries attributes, which are not read")
    items = c.findall("DataItem")
    need = 2 if kind == "side" else 1
    if len(items) < need:
        raise ReadError(
            f"XDMF: set '{name}' of SetType '{set_type}' needs {need} DataItem(s)"
        )

    def ids(item):
        if item.attrib.get("Dimensions") == "0":
            return np.empty(0, dtype=np.int64)
        return np.asarray(reader._read_data_item(item), dtype=np.int64).reshape(-1)

    first = ids(items[0])
    if kind == "side":
        local = ids(items[1]) if len(first) else np.empty(0, dtype=np.int64)
        if len(local) != len(first):
            raise ReadError(
                f"XDMF: set '{name}' has {len(first)} cells but {len(local)} "
                "local face/edge indices"
            )
        entries = np.stack([first, local], axis=1)
    else:
        entries = first
    key = (name, kind)
    if key in regions:
        prev = regions[key]
        entries = np.concatenate([np.asarray(prev.entries), entries])
        dim, tag = prev.dim, prev.tag
    regions[key] = Region(name, kind, entries, dim, tag)


class XdmfReader:
    def __init__(self, filename):
        self.filename = filename

    def read(self):
        parser = ET.XMLParser()
        tree = ET.parse(self.filename, parser)
        root = tree.getroot()
        self._root = root

        if root.tag != "Xdmf":
            raise ReadError()

        version = root.attrib["Version"]

        if version.split(".")[0] == "2":
            return self.read_xdmf2(root)

        if version.split(".")[0] != "3":
            raise ReadError(f"Unknown XDMF version {version}.")

        return self.read_xdmf3(root)

    def _read_data_item(self, data_item, root=None):
        import h5py

        root = root if root is not None else getattr(self, "_root", None)
        visited = set()
        while data_item is not None and "Reference" in data_item.attrib:
            if id(data_item) in visited:
                raise ReadError("XDMF: cyclic DataItem reference")
            visited.add(id(data_item))
            reference = data_item.attrib["Reference"]
            xpath = (
                (data_item.text or "") if reference == "XML" else reference
            ).strip()
            if not xpath.startswith("/") or root is None:
                raise ReadError("XDMF: DataItem reference must be an absolute XPath")
            try:
                targets = root.findall("/".join(xpath.split("/")[2:]))
            except (SyntaxError, KeyError) as exc:
                raise ReadError(f"XDMF: invalid reference XPath {xpath!r}") from exc
            if len(targets) != 1 or targets[0].tag != "DataItem":
                raise ReadError("XDMF: reference must select exactly one DataItem")
            data_item = targets[0]

        if data_item is None:
            raise ReadError("XDMF: missing DataItem")

        dims = [int(d) for d in data_item.get("Dimensions").split()]

        # Actually, `NumberType` is XDMF2 and `DataType` XDMF3, but many files
        # out there use both keys interchangeably.
        if "DataType" in data_item.attrib and "NumberType" in data_item.attrib:
            raise ReadError()
        elif "DataType" in data_item.attrib:
            data_type = data_item.attrib["DataType"]
        elif "NumberType" in data_item.attrib:
            data_type = data_item.attrib["NumberType"]
        else:
            # Default, see
            # <https://xdmf.org/index.php/XDMF_Model_and_Format#XML_Element_.28Xdmf_ClassName.29_and_Default_XML_Attributes>
            data_type = "Float"

        precision = data_item.get("Precision", default="4")
        fmt = data_item.attrib["Format"]

        if fmt == "XML":
            dtype = xdmf_to_numpy_type[(data_type, precision)]
            if data_item.text.strip() == "":
                # https://github.com/numpy/numpy/issues/18435
                data = np.empty((0,), dtype=dtype)
            else:
                data = np.fromstring(data_item.text, dtype=dtype, sep=" ")
            return data.reshape(dims)

        elif fmt == "Binary":
            return np.fromfile(
                data_item.text.strip(), dtype=xdmf_to_numpy_type[(data_type, precision)]
            ).reshape(dims)

        if fmt != "HDF":
            raise ReadError(f"Unknown XDMF Format '{fmt}'.")

        info = data_item.text.strip()
        filename, h5path = info.split(":")

        # The HDF5 file path is given with respect to the XDMF (XML) file.
        dirname = pathlib.Path(self.filename).resolve().parent
        full_hdf5_path = dirname / filename

        f = h5py.File(full_hdf5_path, "r")

        # Some files don't contain the leading slash /.
        if h5path[0] == "/":
            h5path = h5path[1:]

        for key in h5path.split("/"):
            f = f[key]
        # The file's own extent is not bounded by the XML: refuse a dataset
        # larger than its DataItem declares before reading it (the C++ core
        # does the same), rather than allocate whatever a corrupt file claims.
        if dims and int(np.prod(f.shape, dtype=object)) > int(
            np.prod(dims, dtype=object)
        ):
            raise ReadError(
                f"XDMF: HDF dataset '{h5path}' is larger than the Dimensions "
                "its DataItem declares"
            )
        # `[()]` gives a np.ndarray
        return f[()]

    def read_information(self, c_data):
        field_data = {}
        root = ET.fromstring(c_data)
        for child in root:
            str_tag = child.attrib["key"]
            dim = int(child.attrib["dim"])
            if child.text is None:
                raise ReadError()
            num_tag = int(child.text)
            field_data[str_tag] = np.array([num_tag, dim])
        return field_data

    def read_xdmf2(self, root):  # noqa: C901
        domains = list(root)
        if len(domains) != 1:
            raise ReadError()
        domain = domains[0]
        if domain.tag != "Domain":
            raise ReadError()

        grids = list(domain)
        if len(grids) != 1:
            raise ReadError("XDMF reader: Only supports one grid right now.")
        grid = grids[0]
        if grid.tag != "Grid":
            raise ReadError()

        if "GridType" in grid.attrib and grid.attrib["GridType"] != "Uniform":
            raise ReadError()

        points = None
        cells = []
        point_data = {}
        cell_data_raw = {}
        field_data = {}
        regions = {}

        for c in grid:
            if c.tag == "Topology":
                data_items = list(c)
                if len(data_items) != 1:
                    message = (
                        "Need exactly 1 data item in <Topology>, "
                        f"found {len(data_items)}."
                    )
                    if len(data_items) == 0:
                        message += (
                            "\nStructured meshes are not supported, see "
                            "<https://github.com/nschloe/meshio/issues/404>."
                        )
                    raise ReadError(message)

                topology_type = c.attrib["TopologyType"]
                if topology_type == "Mixed":
                    cells = translate_mixed_cells(
                        np.fromstring(
                            data_items[0].text,
                            int,
                            int(data_items[0].get("Dimensions")),
                            " ",
                        )
                    )
                else:
                    data = self._read_data_item(data_items[0], root)
                    cells.append(CellBlock(xdmf_to_meshio_type[topology_type], data))

            elif c.tag == "Geometry":
                if "GeometryType" in c.attrib:
                    geo_type = c.attrib["GeometryType"]
                    if geo_type != "XYZ":
                        raise ReadError(f'Expected GeometryType "XYZ", not {geo_type}.')
                data_items = list(c)
                if len(data_items) != 1:
                    raise ReadError()
                points = self._read_data_item(data_items[0], root)

            elif c.tag == "Information":
                c_data = c.text
                if not c_data:
                    raise ReadError()
                field_data = self.read_information(c_data)

            elif c.tag == "Attribute":
                # assert c.attrib['Active'] == '1'
                # assert c.attrib['AttributeType'] == 'None'

                data_items = list(c)
                if len(data_items) != 1:
                    raise ReadError()

                data = self._read_data_item(data_items[0], root)

                name = c.attrib["Name"]
                if c.attrib["Center"] == "Node":
                    point_data[name] = data
                elif c.attrib["Center"] == "Cell":
                    cell_data_raw[name] = data
                else:
                    # TODO field data?
                    if c.attrib["Center"] != "Grid":
                        raise ReadError()
            elif c.tag == "Set":
                _read_set(self, c, regions)
            else:
                raise ReadError(f"Unknown section '{c.tag}'.")

        cell_data = cell_data_from_raw(cells, cell_data_raw)

        mesh = Mesh(
            points,
            cells,
            point_data=point_data,
            cell_data=cell_data,
            field_data=field_data,
        )
        if regions:
            mesh.regions = list(regions.values())
        return mesh

    def read_xdmf3(self, root):  # noqa: C901
        domains = list(root)
        if len(domains) != 1:
            raise ReadError()
        domain = domains[0]
        if domain.tag != "Domain":
            raise ReadError()

        grids = list(domain)
        if len(grids) != 1:
            raise ReadError("XDMF reader: Only supports one grid right now.")
        grid = grids[0]
        if grid.tag != "Grid":
            raise ReadError()

        points = None
        cells = []
        point_data = {}
        cell_data_raw = {}
        field_data = {}
        regions = {}

        for c in grid:
            if c.tag == "Topology":
                data_items = list(c)
                if len(data_items) != 1:
                    raise ReadError()
                data_item = data_items[0]

                data = self._read_data_item(data_item)

                # The XDMF2 key is `TopologyType`, just `Type` for XDMF3. Allow either
                # one.
                if "Type" in c.attrib and "TopologyType" in c.attrib:
                    raise ReadError()
                elif "Type" in c.attrib:
                    cell_type = c.attrib["Type"]
                else:
                    cell_type = c.attrib["TopologyType"]

                if cell_type == "Mixed":
                    cells = translate_mixed_cells(data)
                else:
                    cells.append(CellBlock(xdmf_to_meshio_type[cell_type], data))

            elif c.tag == "Geometry":
                if "Type" in c.attrib and "GeometryType" in c.attrib:
                    raise ReadError()
                elif "Type" in c.attrib:
                    geometry_type = c.attrib["Type"]
                else:
                    geometry_type = c.attrib["GeometryType"]

                if geometry_type not in ["XY", "XYZ"]:
                    raise ReadError(f'Illegal geometry type "{geometry_type}".')

                data_items = list(c)
                if len(data_items) != 1:
                    raise ReadError()
                data_item = data_items[0]
                points = self._read_data_item(data_item)

            elif c.tag == "Information":
                c_data = c.text
                if not c_data:
                    raise ReadError()
                field_data = self.read_information(c_data)

            elif c.tag == "Attribute":
                # Don't be too strict here: FEniCS, for example, calls this
                # 'AttributeType'.
                # assert c.attrib['Type'] == 'None'

                data_items = list(c)
                if len(data_items) != 1:
                    raise ReadError()
                data_item = data_items[0]

                data = self._read_data_item(data_item)

                if c.attrib["Center"] not in ["Node", "Cell"]:
                    raise ReadError(f"Unknown center '{c.attrib['Center']}'.")

                name = c.attrib["Name"]
                if c.attrib["Center"] == "Node":
                    point_data[name] = data
                else:
                    cell_data_raw[name] = data

            elif c.tag == "Set":
                _read_set(self, c, regions)
            else:
                raise ReadError(f"Unknown section '{c.tag}'.")

        cell_data = cell_data_from_raw(cells, cell_data_raw)

        mesh = Mesh(
            points,
            cells,
            point_data=point_data,
            cell_data=cell_data,
            field_data=field_data,
        )
        if regions:
            mesh.regions = list(regions.values())
        return mesh


class XdmfWriter:
    def __init__(
        self, filename, mesh, data_format="HDF", compression="gzip", compression_opts=4
    ):
        import h5py

        if data_format not in ["XML", "Binary", "HDF"]:
            raise WriteError(
                "Unknown XDMF data format "
                f"'{data_format}' (use 'XML', 'Binary', or 'HDF'.)"
            )

        self.filename = pathlib.Path(filename)
        self.data_format = data_format
        self.data_counter = 0
        self.compression = compression
        self.compression_opts = None if compression is None else compression_opts

        if data_format == "HDF":
            self.h5_filename = self.filename.with_suffix(".h5")
            self.h5_file = h5py.File(self.h5_filename, "w")

        xdmf_file = ET.Element("Xdmf", Version="3.0")

        domain = ET.SubElement(xdmf_file, "Domain")
        grid = ET.SubElement(domain, "Grid", Name="Grid")
        # information = ET.SubElement(
        #     grid, "Information", Name="Information", Value=str(len(mesh.field_data))
        # )

        self.write_points(grid, mesh.points)
        # self.field_data(mesh.field_data, information)
        self.write_cells(mesh.cells, grid)
        self.write_point_data(mesh.point_data, grid)
        self.write_cell_data(mesh.cell_data, grid)
        self.write_sets(mesh, grid)

        ET.register_namespace("xi", "https://www.w3.org/2001/XInclude/")

        write_xml(filename, xdmf_file)

    def numpy_to_xml_string(self, data):
        if self.data_format == "XML":
            s = BytesIO()
            fmt = dtype_to_format_string[data.dtype.name]
            np.savetxt(s, data, fmt)
            return "\n" + s.getvalue().decode()
        elif self.data_format == "Binary":
            base = os.path.splitext(self.filename)[0]
            bin_filename = f"{base}{self.data_counter}.bin"
            self.data_counter += 1
            # write binary data to file
            with open(bin_filename, "wb") as f:
                data.tofile(f)
            return bin_filename

        if self.data_format != "HDF":
            raise WriteError(f'Unknown data format "{self.data_format}"')
        name = f"data{self.data_counter}"
        self.data_counter += 1
        self.h5_file.create_dataset(
            name,
            data=data,
            compression=self.compression,
            compression_opts=self.compression_opts,
        )
        return os.path.basename(self.h5_filename) + ":/" + name

    def _write_ids(self, parent, ids):
        ids = np.asarray(ids, dtype=np.int64).reshape(-1)
        if ids.size == 0:
            ET.SubElement(
                parent,
                "DataItem",
                DataType="Int",
                Dimensions="0",
                Format="XML",
                Precision="8",
            )
            return
        dt, prec = numpy_to_xdmf_dtype[ids.dtype.name]
        data_item = ET.SubElement(
            parent,
            "DataItem",
            DataType=dt,
            Dimensions=str(ids.size),
            Format=self.data_format,
            Precision=prec,
        )
        data_item.text = self.numpy_to_xml_string(ids)

    def _write_set(self, grid, region, set_type, ids, local=None):
        s = ET.SubElement(grid, "Set", Name=region.name, SetType=set_type)
        if region.dim != -1:
            ET.SubElement(s, "Information", Name="meshio++:dim", Value=str(region.dim))
        if region.tag != -1:
            ET.SubElement(s, "Information", Name="meshio++:tag", Value=str(region.tag))
        self._write_ids(s, ids)
        if local is not None:
            self._write_ids(s, local)

    def write_sets(self, mesh, grid):
        """Regions as ``<Set>`` s, the twin of the core's ``xdmf_write_sets``
        (same order: the regions' canonical ``(kind, name, dim, tag)`` order)."""
        order = {"point": 0, "cell": 1, "side": 2}
        regions = sorted(
            getattr(mesh, "regions", []),
            key=lambda r: (order[r.kind], r.name, r.dim, r.tag),
        )
        bases = block_bases(mesh.cells)
        for region in regions:
            if region.kind != "side":
                self._write_set(
                    grid,
                    region,
                    "Node" if region.kind == "point" else "Cell",
                    region.entries,
                )
                continue
            pairs = np.asarray(region.entries, dtype=np.int64).reshape(-1, 2)
            is_edge = np.zeros(len(pairs), dtype=bool)
            for i, cell in enumerate(pairs[:, 0].tolist()):
                b = int(np.searchsorted(bases, cell, side="right")) - 1
                if 0 <= b < len(mesh.cells):
                    is_edge[i] = _block_dim(mesh.cells[b]) == 2
            faces, edges = pairs[~is_edge], pairs[is_edge]
            if len(faces) or not len(edges):
                self._write_set(grid, region, "Face", faces[:, 0], faces[:, 1])
            if len(edges):
                self._write_set(grid, region, "Edge", edges[:, 0], edges[:, 1])

    def write_points(self, grid, points):
        if points.shape[1] > 3:
            raise WriteError("Can only write points up to dimension 3.")

        geometry_type = "XYZ"[: points.shape[1]]

        geo = ET.SubElement(grid, "Geometry", GeometryType=geometry_type)
        dt, prec = numpy_to_xdmf_dtype[points.dtype.name]
        dim = "{} {}".format(*points.shape)
        data_item = ET.SubElement(
            geo,
            "DataItem",
            DataType=dt,
            Dimensions=dim,
            Format=self.data_format,
            Precision=prec,
        )
        data_item.text = self.numpy_to_xml_string(points)

    def write_cells(self, cells, grid):
        if len(cells) == 0:
            return
        if len(cells) == 1:
            meshio_type = cells[0].type
            num_cells = len(cells[0].data)
            xdmf_type = meshio_to_xdmf_type[meshio_type][0]
            topo = ET.SubElement(
                grid,
                "Topology",
                TopologyType=xdmf_type,
                NumberOfElements=str(num_cells),
                NodesPerElement=str(cells[0].data.shape[1]),
            )
            dt, prec = numpy_to_xdmf_dtype[cells[0].data.dtype.name]
            dim = "{} {}".format(*cells[0].data.shape)
            data_item = ET.SubElement(
                topo,
                "DataItem",
                DataType=dt,
                Dimensions=dim,
                Format=self.data_format,
                Precision=prec,
            )
            data_item.text = self.numpy_to_xml_string(cells[0].data)

        else:
            assert len(cells) > 1
            total_num_cells = sum(c.data.shape[0] for c in cells)
            topo = ET.SubElement(
                grid,
                "Topology",
                TopologyType="Mixed",
                NumberOfElements=str(total_num_cells),
            )
            total_num_cell_items = sum(np.prod(c.data.shape) for c in cells)
            num_vertices_and_lines = sum(
                c.data.shape[0] for c in cells if c.type in {"vertex", "line"}
            )
            dim = str(total_num_cell_items + total_num_cells + num_vertices_and_lines)
            cd = np.concatenate(
                [
                    np.hstack(
                        [
                            np.full(
                                (
                                    cell_block.data.shape[0],
                                    2 if cell_block.type in {"vertex", "line"} else 1,
                                ),
                                meshio_type_to_xdmf_index[cell_block.type],
                            ),
                            cell_block.data,
                        ]
                    ).flatten()
                    for cell_block in cells
                ]
            )
            dt, prec = numpy_to_xdmf_dtype[cd.dtype.name]
            data_item = ET.SubElement(
                topo,
                "DataItem",
                DataType=dt,
                Dimensions=dim,
                Format=self.data_format,
                Precision=prec,
            )
            data_item.text = self.numpy_to_xml_string(cd)

    def write_point_data(self, point_data, grid):
        for name, data in point_data.items():
            att = ET.SubElement(
                grid,
                "Attribute",
                Name=name,
                AttributeType=attribute_type(data),
                Center="Node",
            )
            dt, prec = numpy_to_xdmf_dtype[data.dtype.name]
            dim = " ".join([str(s) for s in data.shape])
            data_item = ET.SubElement(
                att,
                "DataItem",
                DataType=dt,
                Dimensions=dim,
                Format=self.data_format,
                Precision=prec,
            )
            data_item.text = self.numpy_to_xml_string(data)

    def write_cell_data(self, cell_data, grid):
        raw = raw_from_cell_data(cell_data)
        for name, data in raw.items():
            att = ET.SubElement(
                grid,
                "Attribute",
                Name=name,
                AttributeType=attribute_type(data),
                Center="Cell",
            )
            dt, prec = numpy_to_xdmf_dtype[data.dtype.name]
            dim = " ".join([str(s) for s in data.shape])
            data_item = ET.SubElement(
                att,
                "DataItem",
                DataType=dt,
                Dimensions=dim,
                Format=self.data_format,
                Precision=prec,
            )
            data_item.text = self.numpy_to_xml_string(data)

    # The original idea was to implement field data as XML CDATA. Unfortunately, in
    # Python's XML, CDATA handled poorly. There are workarounds, e.g.,
    # <https://stackoverflow.com/a/30019607/353337>, but those mess around with
    # ET._serialize_xml and lead to bugs elsewhere
    # <https://github.com/nschloe/meshio/issues/806>.
    # def field_data(self, field_data, information):
    #     info = ET.Element("main")
    #     for name, data in field_data.items():
    #         data_item = ET.SubElement(info, "map", key=name, dim=str(data[1]))
    #         data_item.text = str(data[0])
    #     information.text = ET.CDATA(ET.tostring(info))
    #     information.append(CDATA(ET.tostring(info).decode()))


def write(*args, **kwargs):
    XdmfWriter(*args, **kwargs)
