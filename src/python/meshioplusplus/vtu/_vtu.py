"""
I/O for VTU.
<https://vtk.org/Wiki/VTK_XML_Formats>
<https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf>
"""

import base64
import re
import sys
import zlib

import numpy as np

from .._common import info, join_strings, raw_from_cell_data, replace_space, warn
from .._exceptions import CorruptionError, ReadError
from .._mesh import CellBlock, Mesh
from .._provenance import TAG as _PROVENANCE_TAG
from .._vtk_common import meshio_to_vtk_order, meshio_to_vtk_type, vtk_cells_from_data


class _OptionalCodec:
    """A zlib-like ``compress``/``decompress`` pair backed by an optional module.

    lz4 and zstd are not in the standard library, so the pure-Python reference
    path can only handle them when the matching extra is installed. Failing here
    with a named, actionable error beats an ImportError from deep inside the
    parser -- and beats silently producing a different codec.

    Install with: ``pip install "meshioplusplus[codecs]"``.
    """

    def __init__(self, name: str, extra: str):
        self._name = name
        self._extra = extra

    def _missing(self):
        return ReadError(
            f"{self._name}-compressed VTK XML needs the '{self._extra}' package "
            f'(pip install "meshioplusplus[codecs]"), or a C++ core built with '
            f"-DMESHIOPLUSPLUS_WITH_{self._name.upper()}=ON"
        )

    def compress(self, data, *args, **kwargs):
        raise self._missing()

    def decompress(self, data, *args, **kwargs):
        raise self._missing()


class _Lz4Codec(_OptionalCodec):
    """LZ4 raw *block* format, matching ``vtkLZ4DataCompressor``."""

    def __init__(self):
        super().__init__("lz4", "lz4")

    def compress(self, data, *args, **kwargs):
        try:
            import lz4.block
        except ImportError:
            raise self._missing() from None
        return lz4.block.compress(data, mode="default", store_size=False)

    def decompress(self, data, uncompressed_size=None, **kwargs):
        try:
            import lz4.block
        except ImportError:
            raise self._missing() from None
        if uncompressed_size is None:
            raise ReadError("lz4 block decompression needs the decompressed size")
        return lz4.block.decompress(data, uncompressed_size=uncompressed_size)


class _ZstdCodec(_OptionalCodec):
    """Raw zstd frames, one per VTU block."""

    def __init__(self):
        super().__init__("zstd", "zstandard")

    def compress(self, data, *args, **kwargs):
        try:
            import zstandard
        except ImportError:
            raise self._missing() from None
        return zstandard.ZstdCompressor().compress(data)

    def decompress(self, data, uncompressed_size=None, **kwargs):
        try:
            import zstandard
        except ImportError:
            raise self._missing() from None
        return zstandard.ZstdDecompressor().decompress(data)


_lz4 = _Lz4Codec()
_zstd = _ZstdCodec()

#: ``compression=`` keyword value -> the ``compressor=`` attribute it writes.
_COMPRESSION_TO_ATTR = {
    "zlib": "vtkZLibDataCompressor",
    "lzma": "vtkLZMADataCompressor",
    "lz4": "vtkLZ4DataCompressor",
    "zstd": "vtkZSTDDataCompressor",
}


class _SizeTolerantCodec:
    """Wraps a stdlib codec so it accepts (and ignores) ``uncompressed_size``.

    zlib and lzma record the decompressed size themselves; LZ4's raw block
    format does not. Passing the size unconditionally at the call site keeps
    that difference out of the reader loop.
    """

    def __init__(self, module):
        self._module = module

    def compress(self, data, *args, **kwargs):
        return self._module.compress(data, *args, **kwargs)

    def decompress(self, data, uncompressed_size=None, **kwargs):
        return self._module.decompress(data, **kwargs)


def _compressor_for(name: str):
    """``compressor=`` attribute -> a zlib-like compress/decompress pair.

    Resolved lazily rather than as a module-level dict: ``lzma`` is imported on
    demand elsewhere in this file (it can be missing from a stripped CPython),
    so binding it at import time would break the module entirely.

    lz4 is a real VTK compressor; zstd is a meshio++ extension, since VTK ships
    no ZSTD compressor.
    """
    if name == "vtkZLibDataCompressor":
        return _SizeTolerantCodec(zlib)
    if name == "vtkLZ4DataCompressor":
        return _lz4
    if name == "vtkZSTDDataCompressor":
        return _zstd
    if name == "vtkLZMADataCompressor":
        import lzma

        return _SizeTolerantCodec(lzma)
    raise ReadError(f"Unknown VTK XML compressor '{name}'")


# Paraview 5.8.1's built-in Python doesn't have lzma.
try:
    import lzma
except ModuleNotFoundError:
    lzma = None


def num_bytes_to_num_base64_chars(num_bytes):
    # Rounding up in integer division works by double negation since Python
    # always rounds down.
    return -(-num_bytes // 3) * 4


def _polyhedron_cells_from_data(offsets, faces, faceoffsets, cell_data_raw):
    # In general the number of faces will vary between cells, and the
    # number of nodes vary between faces for each cell. The information
    # will be stored as a List (one item per cell) of lists (one item
    # per face of the cell) of np-arrays of node indices.

    cells = {}
    cell_data = {}

    # The data format for face-cells is:
    # num_faces_cell_0,
    #   num_nodes_face_0, node_ind_0, node_ind_1, ..
    #   num_nodes_face_1, node_ind_0, node_ind_1, ..
    #   ...
    # num_faces_cell_1,
    #   ...
    # See https://vtk.org/Wiki/VTK/Polyhedron_Support for more.

    # The faceoffsets describes the end of the face description for each
    # cell. Switch faceoffsets to give start points, not end points
    faceoffsets = np.append([0], faceoffsets[:-1])

    # Double loop over cells then faces.
    # This will be slow, but seems necessary to cover all cases
    for cell_start in faceoffsets:
        num_faces_this_cell = faces[cell_start]
        faces_this_cell = []
        next_face = cell_start + 1
        for _ in range(num_faces_this_cell):
            num_nodes_this_face = faces[next_face]
            faces_this_cell.append(
                np.array(
                    faces[next_face + 1 : (next_face + num_nodes_this_face + 1)],
                    dtype=int,
                )
            )
            # Increase by number of nodes just read, plus the item giving
            # number of nodes per face
            next_face += num_nodes_this_face + 1

        # Done with this cell
        # Find number of nodes for this cell
        num_nodes_this_cell = np.unique(np.hstack([v for v in faces_this_cell])).size

        key = f"polyhedron{num_nodes_this_cell}"
        if key not in cells.keys():
            cells[key] = []
        cells[key].append(faces_this_cell)

    # The cells will be assigned to blocks according to their number of nodes.
    # This is potentially a reordering, compared to the ordering in faces.
    # Cell data must be reorganized accordingly.

    # Start of the cell-node relations
    start_cn = np.hstack((0, offsets))
    size = np.diff(start_cn)

    # Loop over all cell sizes, find all cells with this size, and store
    # cell data.
    for sz in np.unique(size):
        # Cells with this number of nodes.
        items = np.where(size == sz)[0]

        # Store cell data for this set of cells
        for name, d in cell_data_raw.items():
            if name not in cell_data:
                cell_data[name] = []
            cell_data[name].append(d[items])

    return cells, cell_data


def _organize_cells(point_offsets, cells, cell_data_raw):
    if len(point_offsets) != len(cells):
        raise ReadError("Inconsistent data!")

    out_cells = []

    # IMPLEMENTATION NOTE: The treatment of polyhedral cells is quite a bit different
    # from the other cells; moreover, there are some strong (?) assumptions on such
    # cells. The processing of such cells is therefore moved to a dedicated function for
    # the time being, while all other cell types are treated by the same function.
    # There are still similarities between processing of polyhedral and the rest, so it
    # may be possible to unify the implementations at a later stage.

    # Check if polyhedral cells are present.
    polyhedral_mesh = False
    for c in cells:
        if np.any(c["types"] == 42):  # vtk type 42 is polyhedral
            polyhedral_mesh = True
            break

    if polyhedral_mesh:
        # A single <Piece> only -- multi-piece polyhedral files remain
        # unsupported here, as before.
        if len(cells) > 1:
            raise ValueError("Implementation assumes single set of cells")

        types = cells[0]["types"].ravel()
        if np.all(types == 42):
            # Every cell is a polyhedron: the historical fast path, unchanged.
            cls, cell_data = _polyhedron_cells_from_data(
                cells[0]["offsets"],
                cells[0]["faces"],
                cells[0]["faceoffsets"],
                cell_data_raw[0],
            )
            for tp, c in cls.items():
                out_cells.append(CellBlock(tp, c))
        else:
            # MIXED. `faceoffsets` carries -1 for a cell that is not a
            # polyhedron, which is exactly how VTU expresses this -- an
            # OpenFOAM mesh always mixes hexahedra, polyhedra and boundary
            # faces. Before v9.19.0 this raised "Cannot handle combinations of
            # polyhedra with other cells"; the C++ reader handles it, and this
            # twin must agree or a no-core build (Windows CI) would refuse a
            # file every other build reads.
            offsets = cells[0]["offsets"].ravel()
            conn = cells[0]["connectivity"].ravel()
            faces = np.asarray(cells[0]["faces"]).ravel()
            faceoffsets = np.asarray(cells[0]["faceoffsets"]).ravel()
            cell_data = {}
            # Walk contiguous runs of same-"is it a polyhedron" so the output
            # block order follows the file's own cell order.
            is_poly = types == 42
            start = 0
            n = len(types)
            while start < n:
                end = start + 1
                while end < n and is_poly[end] == is_poly[start]:
                    end += 1
                sub_raw = {k: v[start:end] for k, v in cell_data_raw[0].items()}
                if is_poly[start]:
                    cls, cd = _polyhedron_cells_from_data(
                        offsets[start:end],
                        faces,
                        faceoffsets[start:end],
                        sub_raw,
                    )
                    for tp, c in cls.items():
                        out_cells.append(CellBlock(tp, c))
                else:
                    first_node = 0 if start == 0 else offsets[start - 1]
                    blocks, cd = vtk_cells_from_data(
                        conn[first_node : offsets[end - 1]],
                        offsets[start:end] - first_node,
                        types[start:end],
                        sub_raw,
                    )
                    for c in blocks:
                        out_cells.append(CellBlock(c.type, c.data + first_node))
                for k, v in cd.items():
                    cell_data.setdefault(k, []).extend(v)
                start = end

    else:
        for offset, cls, cdr in zip(point_offsets, cells, cell_data_raw):
            cls, cell_data = vtk_cells_from_data(
                cls["connectivity"].ravel(),
                cls["offsets"].ravel(),
                cls["types"].ravel(),
                cdr,
            )

        for c in cls:
            out_cells.append(CellBlock(c.type, c.data + offset))

    return out_cells, cell_data


def get_grid(root):
    grid = None
    appended_data = None
    for c in root:
        if c.tag == "UnstructuredGrid":
            if grid is not None:
                raise ReadError("More than one UnstructuredGrid found.")
            grid = c
        else:
            if c.tag != "AppendedData":
                raise ReadError(f"Unknown main tag '{c.tag}'.")
            if appended_data is not None:
                raise ReadError("More than one AppendedData section found.")
            if c.attrib["encoding"] != "base64":
                raise ReadError("")
            appended_data = c.text.strip()
            # The appended data always begins with a (meaningless) underscore.
            if appended_data[0] != "_":
                raise ReadError()
            appended_data = appended_data[1:]

    if grid is None:
        raise ReadError("No UnstructuredGrid found.")
    return grid, appended_data


def _parse_raw_binary(filename):
    from xml.etree import ElementTree as ET

    with open(filename, "rb") as f:
        raw = f.read()

    try:
        res = re.search(re.compile(b'<AppendedData[^>]+(?:">)'), raw)
        assert res is not None
        i_start = res.end()
        i_stop = raw.find(b"</AppendedData>")
    except Exception:
        raise ReadError()

    header = raw[:i_start].decode()
    footer = raw[i_stop:].decode()
    data = raw[i_start:i_stop].split(b"_", 1)[1].rsplit(b"\n", 1)[0]

    root = ET.fromstring(header + footer)

    dtype = vtu_to_numpy_type[root.get("header_type", "UInt32")]
    if "byte_order" in root.attrib:
        dtype = dtype.newbyteorder(
            "<" if root.get("byte_order") == "LittleEndian" else ">"
        )

    appended_data_tag = root.find("AppendedData")
    assert appended_data_tag is not None
    appended_data_tag.set("encoding", "base64")

    compressor = root.get("compressor")
    if compressor is None:
        arrays = ""
        i = 0
        while i < len(data):
            # The following find() runs into issues if offset is padded with spaces, see
            # <https://github.com/nschloe/meshio/issues/1135>. It works in ParaView.
            # Unfortunately, Python's built-in XML tree can't handle regexes, see
            # <https://stackoverflow.com/a/38810731/353337>.
            da_tag = root.find(f".//DataArray[@offset='{i}']")
            if da_tag is None:
                raise RuntimeError(f"Could not find .//DataArray[@offset='{i}']")
            da_tag.set("offset", str(len(arrays)))

            block_size = int(np.frombuffer(data[i : i + dtype.itemsize], dtype)[0])
            arrays += base64.b64encode(
                data[i : i + block_size + dtype.itemsize]
            ).decode()
            i += block_size + dtype.itemsize

    else:
        c = _compressor_for(compressor)
        root.attrib.pop("compressor")

        # raise ReadError("Compressed raw binary VTU files not supported.")
        arrays = ""
        i = 0
        while i < len(data):
            da_tag = root.find(f".//DataArray[@offset='{i}']")
            assert da_tag is not None
            da_tag.set("offset", str(len(arrays)))

            num_blocks = int(np.frombuffer(data[i : i + dtype.itemsize], dtype)[0])
            num_header_items = 3 + num_blocks
            num_header_bytes = num_header_items * dtype.itemsize
            header = np.frombuffer(data[i : i + num_header_bytes], dtype)

            block_data = b""
            j = 0
            for k in range(num_blocks):
                block_size = int(header[k + 3])
                block_data += c.decompress(
                    data[
                        i + j + num_header_bytes : i + j + block_size + num_header_bytes
                    ]
                )
                j += block_size

            block_size = np.array([len(block_data)]).astype(dtype).tobytes()
            arrays += base64.b64encode(block_size + block_data).decode()

            i += j + num_header_bytes

    appended_data_tag.text = "_" + arrays
    return root


vtu_to_numpy_type = {
    "Float32": np.dtype(np.float32),
    "Float64": np.dtype(np.float64),
    "Int8": np.dtype(np.int8),
    "Int16": np.dtype(np.int16),
    "Int32": np.dtype(np.int32),
    "Int64": np.dtype(np.int64),
    "UInt8": np.dtype(np.uint8),
    "UInt16": np.dtype(np.uint16),
    "UInt32": np.dtype(np.uint32),
    "UInt64": np.dtype(np.uint64),
}
numpy_to_vtu_type = {v: k for k, v in vtu_to_numpy_type.items()}


class VtuReader:
    """Helper class for reading VTU files. Some properties are global to the file (e.g.,
    byte_order), and instead of passing around these parameters, make them properties of
    this class.
    """

    def __init__(self, filename):  # noqa: C901
        from xml.etree import ElementTree as ET

        parser = ET.XMLParser()
        try:
            tree = ET.parse(str(filename), parser)
            root = tree.getroot()
        except ET.ParseError:
            root = _parse_raw_binary(str(filename))

        if root.tag != "VTKFile":
            raise ReadError(f"Expected tag 'VTKFile', found {root.tag}")
        if root.attrib["type"] != "UnstructuredGrid":
            tpe = root.attrib["type"]
            raise ReadError(f"Expected type UnstructuredGrid, found {tpe}")

        if "version" in root.attrib:
            version = root.attrib["version"]
            if version not in ["0.1", "1.0"]:
                raise ReadError(f"Unknown VTU file version '{version}'.")

        # fix empty NumberOfComponents attributes as produced by Firedrake
        for da_tag in root.findall(".//DataArray[@NumberOfComponents='']"):
            da_tag.attrib.pop("NumberOfComponents")

        if "compressor" in root.attrib:
            assert root.attrib["compressor"] in [
                "vtkLZMADataCompressor",
                "vtkZLibDataCompressor",
                "vtkLZ4DataCompressor",
                "vtkZSTDDataCompressor",
            ]
            self.compression = root.attrib["compressor"]
        else:
            self.compression = None

        self.header_type = (
            root.attrib["header_type"] if "header_type" in root.attrib else "UInt32"
        )

        try:
            self.byte_order = root.attrib["byte_order"]
            if self.byte_order not in ["LittleEndian", "BigEndian"]:
                raise ReadError(f"Unknown byte order '{self.byte_order}'.")
        except KeyError:
            self.byte_order = None

        grid, self.appended_data = get_grid(root)

        pieces = []
        field_data = {}
        for c in grid:
            if c.tag == "Piece":
                pieces.append(c)
            elif c.tag == "FieldData":
                # TODO test field data
                for data_array in c:
                    field_data[data_array.attrib["Name"]] = self.read_data(data_array)
            else:
                raise ReadError(f"Unknown grid subtag '{c.tag}'.")

        if not pieces:
            raise ReadError("No Piece found.")

        points = []
        cells = []
        point_data = []
        cell_data_raw = []

        for piece in pieces:
            piece_cells = {}
            piece_point_data = {}
            piece_cell_data_raw = {}

            num_points = int(piece.attrib["NumberOfPoints"])
            num_cells = int(piece.attrib["NumberOfCells"])

            for child in piece:
                if child.tag == "Points":
                    data_arrays = list(child)
                    if len(data_arrays) != 1:
                        raise ReadError()
                    data_array = data_arrays[0]

                    if data_array.tag != "DataArray":
                        raise ReadError()

                    pts = self.read_data(data_array)

                    num_components = int(data_array.attrib["NumberOfComponents"])
                    points.append(pts.reshape(num_points, num_components))

                elif child.tag == "Cells":
                    for data_array in child:
                        if data_array.tag != "DataArray":
                            raise ReadError()
                        piece_cells[data_array.attrib["Name"]] = self.read_data(
                            data_array
                        )

                    if len(piece_cells["offsets"]) != num_cells:
                        raise ReadError()
                    if len(piece_cells["types"]) != num_cells:
                        raise ReadError()

                    cells.append(piece_cells)

                elif child.tag == "PointData":
                    for c in child:
                        if c.tag != "DataArray":
                            raise ReadError()
                        try:
                            piece_point_data[c.attrib["Name"]] = self.read_data(c)
                        except CorruptionError as e:
                            warn(e.args[0] + " Skipping.")

                    point_data.append(piece_point_data)

                elif child.tag == "CellData":
                    for c in child:
                        if c.tag != "DataArray":
                            raise ReadError()
                        piece_cell_data_raw[c.attrib["Name"]] = self.read_data(c)

                    cell_data_raw.append(piece_cell_data_raw)
                else:
                    raise ReadError(f"Unknown tag '{child.tag}'.")

        if not cell_data_raw:
            cell_data_raw = [{}] * len(cells)

        if len(cell_data_raw) != len(cells):
            raise ReadError()

        point_offsets = np.cumsum([0] + [pts.shape[0] for pts in points][:-1])

        # Now merge across pieces
        if not points:
            raise ReadError()
        self.points = np.concatenate(points)

        if point_data:
            self.point_data = {
                key: np.concatenate([pd[key] for pd in point_data])
                for key in point_data[0]
            }
        else:
            self.point_data = None

        self.cells, self.cell_data = _organize_cells(
            point_offsets, cells, cell_data_raw
        )
        self.field_data = field_data

    def read_uncompressed_binary(self, data, dtype):
        byte_string = base64.b64decode(data)

        # the first item is the total_num_bytes, given in header_dtype
        header_dtype = vtu_to_numpy_type[self.header_type]
        if self.byte_order is not None:
            header_dtype = header_dtype.newbyteorder(
                "<" if self.byte_order == "LittleEndian" else ">"
            )
        num_header_bytes = np.dtype(header_dtype).itemsize
        total_num_bytes = np.frombuffer(byte_string[:num_header_bytes], header_dtype)[0]

        # Check if block size was decoded separately
        # (so decoding stopped after block size due to padding)
        if len(byte_string) == num_header_bytes:
            header_len = len(base64.b64encode(byte_string))
            byte_string = base64.b64decode(data[header_len:])
        else:
            byte_string = byte_string[num_header_bytes:]

        # Read the block data; multiple blocks possible here?
        if self.byte_order is not None:
            dtype = dtype.newbyteorder(
                "<" if self.byte_order == "LittleEndian" else ">"
            )
        return np.frombuffer(byte_string[:total_num_bytes], dtype=dtype)

    def read_compressed_binary(self, data, dtype):
        # first read the block size; it determines the size of the header
        header_dtype = vtu_to_numpy_type[self.header_type]
        if self.byte_order is not None:
            header_dtype = header_dtype.newbyteorder(
                "<" if self.byte_order == "LittleEndian" else ">"
            )
        num_bytes_per_item = np.dtype(header_dtype).itemsize
        num_chars = num_bytes_to_num_base64_chars(num_bytes_per_item)
        byte_string = base64.b64decode(data[:num_chars])[:num_bytes_per_item]
        num_blocks = np.frombuffer(byte_string, header_dtype)[0]

        # read the entire header
        num_header_items = 3 + int(num_blocks)
        num_header_bytes = num_bytes_per_item * num_header_items
        num_header_chars = num_bytes_to_num_base64_chars(num_header_bytes)
        byte_string = base64.b64decode(data[:num_header_chars])
        header = np.frombuffer(byte_string, header_dtype)

        # num_blocks = header[0]
        max_uncompressed_block_size = int(header[1])
        last_uncompressed_block_size = int(header[2])
        block_sizes = header[3:]

        # Read the block data
        byte_array = base64.b64decode(data[num_header_chars:])
        if self.byte_order is not None:
            dtype = dtype.newbyteorder(
                "<" if self.byte_order == "LittleEndian" else ">"
            )

        byte_offsets = np.empty(block_sizes.shape[0] + 1, dtype=block_sizes.dtype)
        byte_offsets[0] = 0
        np.cumsum(block_sizes, out=byte_offsets[1:])

        assert self.compression is not None
        c = _compressor_for(self.compression)  # already the compressor= attribute

        # process the compressed data
        # Every block decompresses to max_uncompressed_block_size except the
        # last. zlib/lzma infer that themselves, but LZ4's raw block format
        # cannot -- it has no size field -- so it is passed explicitly. The VTU
        # header carries it, which is exactly why the raw block format is
        # usable here at all.
        block_data = np.concatenate(
            [
                np.frombuffer(
                    c.decompress(
                        byte_array[byte_offsets[k] : byte_offsets[k + 1]],
                        uncompressed_size=(
                            last_uncompressed_block_size
                            if k + 1 == num_blocks
                            else max_uncompressed_block_size
                        ),
                    ),
                    dtype=dtype,
                )
                for k in range(num_blocks)
            ]
        )

        return block_data

    def read_data(self, c):
        fmt = c.attrib["format"] if "format" in c.attrib else "ascii"

        data_type = c.attrib["type"]
        try:
            dtype = vtu_to_numpy_type[data_type]
        except KeyError:
            raise ReadError(f"Illegal data type '{data_type}'.")

        if fmt == "ascii":
            # ascii
            if c.text.strip() == "":
                # https://github.com/numpy/numpy/issues/18435
                data = np.empty((0,), dtype=dtype)
            else:
                data = np.fromstring(c.text, dtype=dtype, sep=" ")
        elif fmt == "binary":
            reader = (
                self.read_uncompressed_binary
                if self.compression is None
                else self.read_compressed_binary
            )
            data = reader(c.text.strip(), dtype)
        elif fmt == "appended":
            offset = int(c.attrib["offset"])
            reader = (
                self.read_uncompressed_binary
                if self.compression is None
                else self.read_compressed_binary
            )
            assert self.appended_data is not None
            data = reader(self.appended_data[offset:], dtype)
        else:
            raise ReadError(f"Unknown data format '{fmt}'.")

        if "NumberOfComponents" in c.attrib:
            nc = int(c.attrib["NumberOfComponents"])
            try:
                data = data.reshape(-1, nc)
            except ValueError:
                name = c.attrib["Name"]
                raise CorruptionError(
                    "VTU file corrupt. "
                    + f"The size of the data array '{name}' is {data.size} "
                    + f"which doesn't fit the number of components {nc}."
                )
        return data


def read(filename):
    reader = VtuReader(filename)
    return Mesh(
        reader.points,
        reader.cells,
        point_data=reader.point_data,
        cell_data=reader.cell_data,
        field_data=reader.field_data,
    )


def _chunk_it(array, n):
    k = 0
    while k * n < len(array):
        yield array[k * n : (k + 1) * n]
        k += 1


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    # Writing XML with an etree required first transforming the (potentially large)
    # arrays into string, which are much larger in memory still. This makes this writer
    # very memory hungry. See <https://stackoverflow.com/q/59272477/353337>.
    from .._cxml import etree as ET

    # Check if the mesh contains polyhedral cells, this will require special treatment
    # in certain places.
    is_polyhedron_grid = False
    for c in mesh.cells:
        if c.type.startswith("polyhedron"):
            is_polyhedron_grid = True
            break

    # Mixing polyhedra with other cell types IS expressible: `faceoffsets`
    # carries -1 for a non-polyhedral cell. See the reader above and
    # doc/polyhedra.md; the C++ writer does the same.

    if not binary:
        warn("VTU ASCII files are only meant for debugging.")

    if mesh.points.shape[1] == 2:
        warn(
            "VTU requires 3D points, but 2D points given. "
            "Appending 0 third component."
        )
        points = np.column_stack([mesh.points, np.zeros_like(mesh.points[:, 0])])
    else:
        points = mesh.points

    if mesh.point_sets:
        info(
            "VTU format cannot write point_sets. Converting them to point_data...",
            highlight=False,
        )
        key, _ = join_strings(list(mesh.point_sets.keys()))
        key, _ = replace_space(key)
        mesh.point_sets_to_data(key)

    if mesh.cell_sets:
        info(
            "VTU format cannot write cell_sets. Converting them to cell_data...",
            highlight=False,
        )
        key, _ = join_strings(list(mesh.cell_sets.keys()))
        key, _ = replace_space(key)
        mesh.cell_sets_to_data(key)

    vtk_file = ET.Element(
        "VTKFile",
        type="UnstructuredGrid",
        version="0.1",
        # Use the native endianness. Not strictly necessary, but this simplifies things
        # a bit.
        byte_order=("LittleEndian" if sys.byteorder == "little" else "BigEndian"),
    )

    if header_type is None:
        header_type = "UInt32"
    else:
        vtk_file.set("header_type", header_type)
    assert header_type is not None

    if binary and compression:
        # TODO lz4 <https://vtk.org/doc/nightly/html/classvtkDataCompressor.html>
        compressions = {
            "lzma": "vtkLZMADataCompressor",
            "zlib": "vtkZLibDataCompressor",
            "lz4": "vtkLZ4DataCompressor",
            "zstd": "vtkZSTDDataCompressor",
        }
        assert compression in compressions
        vtk_file.set("compressor", compressions[compression])

    # swap the data to match the system byteorder
    # Don't use byteswap to make sure that the dtype is changed; see
    # <https://github.com/numpy/numpy/issues/10372>.
    points = points.astype(points.dtype.newbyteorder("="), copy=False)
    for k, cell_block in enumerate(mesh.cells):
        cell_type = cell_block.type
        data = cell_block.data
        # Treatment of polyhedra is different from other types
        if is_polyhedron_grid:
            new_cell_info = []
            for cell_info in data:
                new_face_info = []
                for face_info in cell_info:
                    face_info = np.asarray(face_info)
                    new_face_info.append(
                        face_info.astype(face_info.dtype.newbyteorder("="), copy=False)
                    )
                new_cell_info.append(new_face_info)
            mesh.cells[k] = CellBlock(cell_type, new_cell_info)
        else:
            mesh.cells[k] = CellBlock(
                cell_type, data.astype(data.dtype.newbyteorder("="), copy=False)
            )
    for key, data in mesh.point_data.items():
        mesh.point_data[key] = data.astype(data.dtype.newbyteorder("="), copy=False)

    for data in mesh.cell_data.values():
        for k, dat in enumerate(data):
            data[k] = dat.astype(dat.dtype.newbyteorder("="), copy=False)
    for key, data in mesh.field_data.items():
        mesh.field_data[key] = data.astype(data.dtype.newbyteorder("="), copy=False)

    def numpy_to_xml_array(parent, name, data):
        vtu_type = numpy_to_vtu_type[data.dtype]
        fmt = "{:.11e}" if vtu_type.startswith("Float") else "{:d}"
        da = ET.SubElement(parent, "DataArray", type=vtu_type, Name=name)
        if len(data.shape) == 2:
            da.set("NumberOfComponents", f"{data.shape[1]}")

        def text_writer_compressed(f):
            max_block_size = 32768
            data_bytes = data.tobytes()

            # round up
            num_blocks = -int(-len(data_bytes) // max_block_size)
            last_block_size = len(data_bytes) - (num_blocks - 1) * max_block_size

            # It's too bad that we have to keep all blocks in memory. This is
            # necessary because the header, written first, needs to know the
            # lengths of all blocks. Also, the blocks are encoded _after_ having
            # been concatenated.
            c = _compressor_for(_COMPRESSION_TO_ATTR[compression])
            compressed_blocks = [
                # This compress is the slowest part of the writer
                c.compress(block)
                for block in _chunk_it(data_bytes, max_block_size)
            ]

            # collect header
            header = np.array(
                [num_blocks, max_block_size, last_block_size]
                + [len(b) for b in compressed_blocks],
                dtype=vtu_to_numpy_type[header_type],
            )
            f.write(base64.b64encode(header.tobytes()).decode())
            f.write(base64.b64encode(b"".join(compressed_blocks)).decode())

        def text_writer_uncompressed(f):
            data_bytes = data.tobytes()
            # collect header
            header = np.array(len(data_bytes), dtype=vtu_to_numpy_type[header_type])
            f.write(base64.b64encode(header.tobytes() + data_bytes).decode())

        def text_writer_ascii(f):
            # This write() loop is the bottleneck for the write. Alternatives:
            # savetxt is super slow:
            #   np.savetxt(f, data.reshape(-1), fmt=fmt)
            # joining and writing is a bit faster, but consumes huge amounts of
            # memory:
            #   f.write("\n".join(map(fmt.format, data.reshape(-1))))
            for item in data.reshape(-1):
                f.write((fmt + "\n").format(item))

        if binary:
            da.set("format", "binary")
            da.text_writer = (
                text_writer_compressed if compression else text_writer_uncompressed
            )
        else:
            da.set("format", "ascii")
            da.text_writer = text_writer_ascii

    def _polyhedron_face_cells(face_cells):
        # Define the faces of each cell on the format specified for VTU Polyhedron
        # cells. These are defined in Mesh.polyhedron_faces, as block data. The block
        # consists of a nested list (outer list represents cell, inner is faces for this
        # cells), where the items of the inner list are the nodes of specific faces.
        #
        # The output format is specified at https://vtk.org/Wiki/VTK/Polyhedron_Support

        # Initialize array for size of data per cell.
        data_size_per_cell = np.zeros(len(face_cells), dtype=int)

        # The data itself is of unknown size, and cannot be initialized
        data = []
        for ci, cell in enumerate(face_cells):
            # Number of faces for this cell
            data.append(len(cell))
            for face in cell:
                # Number of nodes for this face
                data.append(face.size)
                # The nodes themselves
                data += face.tolist()

            data_size_per_cell[ci] = len(data)

        # The returned data corresponds to the faces and faceoffsets fields in the
        # vtu polyhedron data format
        return data, data_size_per_cell.tolist()

    comment = ET.Comment(_PROVENANCE_TAG)
    vtk_file.insert(1, comment)

    grid = ET.SubElement(vtk_file, "UnstructuredGrid")

    total_num_cells = sum(len(c.data) for c in mesh.cells)
    piece = ET.SubElement(
        grid,
        "Piece",
        NumberOfPoints=f"{len(points)}",
        NumberOfCells=f"{total_num_cells}",
    )

    # points
    if points is not None:
        pts = ET.SubElement(piece, "Points")
        numpy_to_xml_array(pts, "Points", points)

    if mesh.cells is not None and len(mesh.cells) > 0:
        cls = ET.SubElement(piece, "Cells")

        faces = None
        faceoffsets = None

        if is_polyhedron_grid:
            # The VTK polyhedron format requires both Cell-node connectivity, and a
            # definition of faces. The cell-node relation must be recoved from the
            # cell-face-nodes currently in CellBlocks.
            # NOTE: If polyhedral cells are implemented for more mesh types, this code
            # block may be useful for those as well.
            con = []
            num_nodes_per_cell = []
            for block in mesh.cells:
                if not block.type.startswith("polyhedron"):
                    # A non-polyhedral block in a mixed grid keeps its ordinary
                    # connectivity; only its `faceoffsets` entry is special
                    # (-1, written below). Before v9.19.0 this loop assumed
                    # every block was polyhedral, because mixing was refused.
                    d = block.data
                    new_order = meshio_to_vtk_order(block.type)
                    if new_order is not None:
                        d = d[:, new_order]
                    for row in d:
                        con += row.tolist()
                        num_nodes_per_cell.append(len(row))
                    continue
                for cell in block.data:
                    nodes_this_cell = []
                    for face in cell:
                        nodes_this_cell += np.asarray(face).tolist()
                    unique_nodes = np.unique(nodes_this_cell).tolist()

                    con += unique_nodes
                    num_nodes_per_cell.append(len(unique_nodes))

            connectivity = np.array(con)
            # offsets = np.hstack(([0], np.cumsum(num_nodes_per_cell)[:-1]))
            offsets = np.cumsum(num_nodes_per_cell)

            # Initialize data structures for polyhedral cells
            faces = []
            faceoffsets = []

        else:
            # create connectivity, offset, type arrays
            connectivity = []
            for v in mesh.cells:
                d = v.data
                new_order = meshio_to_vtk_order(v.type)
                if new_order is not None:
                    d = d[:, new_order]
                connectivity.append(d.flatten())
            connectivity = np.concatenate(connectivity)

            # offset (points to the first element of the next cell)
            offsets = [
                v.data.shape[1]
                * np.arange(1, v.data.shape[0] + 1, dtype=connectivity.dtype)
                for v in mesh.cells
            ]
            for k in range(1, len(offsets)):
                offsets[k] += offsets[k - 1][-1]
            offsets = np.concatenate(offsets)

        # types
        types_array = []
        for cell_block in mesh.cells:
            key = cell_block.type
            # some adaptions for polyhedron
            if key.startswith("polyhedron"):
                # Get face-cell relation on the vtu format. See comments in helper
                # function for more information of how to specify this.
                faces_loc, faceoffsets_loc = _polyhedron_face_cells(cell_block.data)
                # Adjust offsets to global numbering
                assert faceoffsets is not None
                if len(faceoffsets) > 0:
                    faceoffsets_loc = [fi + faceoffsets[-1] for fi in faceoffsets_loc]

                assert faces is not None
                faces += faces_loc
                faceoffsets += faceoffsets_loc
                key = "polyhedron"
            elif is_polyhedron_grid:
                # -1 per cell: the format's own way of saying "this cell is not
                # a polyhedron, it has no entry in the faces stream". Without
                # these the faceoffsets array would be shorter than the cell
                # count and every polyhedron after the first ordinary block
                # would be mis-located.
                assert faceoffsets is not None
                faceoffsets += [-1] * len(cell_block)

            types_array.append(np.full(len(cell_block), meshio_to_vtk_type[key]))

        types = np.concatenate(
            types_array
            # [np.full(len(v), meshio_to_vtk_type[k]) for k, v in mesh.cells]
        )

        numpy_to_xml_array(cls, "connectivity", connectivity)
        numpy_to_xml_array(cls, "offsets", offsets)
        numpy_to_xml_array(cls, "types", types)

        if is_polyhedron_grid:
            # Also store face-node relation
            numpy_to_xml_array(cls, "faces", np.array(faces, dtype=int))
            numpy_to_xml_array(cls, "faceoffsets", np.array(faceoffsets, dtype=int))

    if mesh.point_data:
        pd = ET.SubElement(piece, "PointData")
        for name, data in mesh.point_data.items():
            numpy_to_xml_array(pd, name, data)

    if mesh.cell_data:
        cd = ET.SubElement(piece, "CellData")
        for name, data in raw_from_cell_data(mesh.cell_data).items():
            numpy_to_xml_array(cd, name, data)

    # write_xml(filename, vtk_file, pretty_xml)
    tree = ET.ElementTree(vtk_file)
    tree.write(filename)


# NOTE: format registration now lives in meshioplusplus/vtu/__init__.py, which wraps the
# reader/writer below with the C++-backed fast paths.
