"""GIFTI surfaces: the Python reference of ``formats/gifti.cpp`` (roadmap 2.4.4). See
that header for what is read and written; both engines write the same bytes."""

from __future__ import annotations

import base64
import xml.etree.ElementTree as ET
import zlib

import numpy as np

from .._common import warn
from .._exceptions import ReadError, WriteError
from .._mesh import Mesh
from .._regions import Region

_TYPES = {
    "NIFTI_TYPE_UINT8": np.uint8,
    "NIFTI_TYPE_INT8": np.int8,
    "NIFTI_TYPE_INT16": np.int16,
    "NIFTI_TYPE_UINT16": np.uint16,
    "NIFTI_TYPE_INT32": np.int32,
    "NIFTI_TYPE_UINT32": np.uint32,
    "NIFTI_TYPE_INT64": np.int64,
    "NIFTI_TYPE_UINT64": np.uint64,
    "NIFTI_TYPE_FLOAT32": np.float32,
    "NIFTI_TYPE_FLOAT64": np.float64,
}


def _meta(node, key):
    md = node.find("MetaData")
    if md is None:
        return ""
    for item in md.findall("MD"):
        if (item.findtext("Name") or "") == key:
            return item.findtext("Value") or ""
    return ""


def _decode(array):
    name = array.get("DataType", "")
    if name not in _TYPES:
        raise ReadError(f"meshio++: gifti: DataType '{name}' is not supported")
    dtype = np.dtype(_TYPES[name])
    try:
        rank = int(array.get("Dimensionality", "1"))
    except ValueError:
        rank = -1
    if not 1 <= rank <= 4:
        raise ReadError("meshio++: gifti: Dimensionality must be 1 to 4")
    shape = []
    for d in range(rank):
        try:
            n = int(array.get(f"Dim{d}", "-1"))
        except ValueError:
            n = -1
        if n < 0 or n > (1 << 40):
            raise ReadError(f"meshio++: gifti: Dim{d} is missing or invalid")
        shape.append(n)
    count = int(np.prod(shape)) if shape else 1
    encoding = array.get("Encoding", "")
    text = array.findtext("Data") or ""
    if encoding == "ASCII":
        try:
            values = np.asarray(text.split()[:count], dtype=np.float64)
        except ValueError:
            raise ReadError(
                "meshio++: gifti: an ASCII DataArray holds a non-number"
            ) from None
        if values.size < count:
            raise ReadError(
                "meshio++: gifti: an ASCII DataArray holds fewer values than its dimensions"
            )
        data = values.astype(dtype)
    elif encoding in ("Base64Binary", "GZipBase64Binary"):
        raw = base64.b64decode("".join(text.split()) + "===", validate=False)
        if encoding == "GZipBase64Binary":
            try:
                raw = zlib.decompress(raw)
            except zlib.error as exc:
                raise ReadError(f"meshio++: gifti: {exc}") from None
        if len(raw) != count * dtype.itemsize:
            raise ReadError(
                f"meshio++: gifti: a DataArray holds {len(raw)} bytes, not the "
                f"{count * dtype.itemsize} its dimensions need"
            )
        order = ">" if array.get("Endian") == "BigEndian" else "<"
        data = np.frombuffer(raw, dtype=dtype.newbyteorder(order)).astype(dtype)
    else:
        raise ReadError(
            f"meshio++: gifti: Encoding '{encoding}' is not read (ExternalFileBinary "
            "needs its sibling file)"
        )
    if array.get("ArrayIndexingOrder") == "ColumnMajorOrder" and rank == 2:
        data = data.reshape(shape[::-1]).T
    return np.ascontiguousarray(data.reshape(shape))


def read(filename):
    try:
        root = ET.parse(str(filename)).getroot()
    except (ET.ParseError, OSError) as exc:
        raise ReadError(
            f"meshio++: gifti: could not parse {filename} as XML: {exc}"
        ) from None
    if root.tag != "GIFTI":
        raise ReadError(f"meshio++: gifti: '{filename}' has no <GIFTI> root")
    points, cells, field_data = None, [], {}
    later = []
    for index, array in enumerate(root.findall("DataArray")):
        intent = array.get("Intent", "")
        data = _decode(array)
        if intent == "NIFTI_INTENT_POINTSET" and points is None:
            if data.ndim != 2 or data.shape[1] != 3:
                raise ReadError("meshio++: gifti: the POINTSET is not (n, 3)")
            points = data.astype(np.float64)
            xform = array.find("CoordinateSystemTransformMatrix")
            if xform is not None:
                try:
                    m = np.asarray(
                        (xform.findtext("MatrixData") or "").split()[:16],
                        dtype=np.float64,
                    )
                except ValueError:
                    m = np.zeros(0)
                if m.size != 16:
                    raise ReadError("meshio++: gifti: MatrixData needs 16 numbers")
                field_data["gifti:transform"] = m.reshape(4, 4)
                for key, tag in (
                    ("gifti:data_space", "DataSpace"),
                    ("gifti:transformed_space", "TransformedSpace"),
                ):
                    field_data[key] = np.frombuffer(
                        (xform.findtext(tag) or "").strip().encode(), dtype=np.uint8
                    ).copy()
        elif intent == "NIFTI_INTENT_TRIANGLE":
            if data.ndim != 2 or data.shape[1] != 3:
                raise ReadError("meshio++: gifti: the TRIANGLE array is not (n, 3)")
            later.append(("\x01triangle", data.astype(np.int64), False))
        else:
            name = _meta(array, "Name") or (
                "gifti:label"
                if intent == "NIFTI_INTENT_LABEL"
                else f"gifti:{intent}_{index}"
            )
            later.append((name, data, intent == "NIFTI_INTENT_LABEL"))
    if points is None:
        raise ReadError(
            f"meshio++: gifti: '{filename}' has no NIFTI_INTENT_POINTSET array"
        )
    np_ = len(points)
    point_data, label_name = {}, None
    for name, data, is_label in later:
        if name == "\x01triangle":
            if data.size and (data.min() < 0 or data.max() >= np_):
                raise ReadError(
                    "meshio++: gifti: a triangle names a point out of range"
                )
            cells.append(("triangle", data))
            continue
        if data.ndim == 0 or data.shape[0] != np_:
            warn(f"gifti: array '{name}' does not hold one row per point; skipped")
            continue
        if data.ndim > 2:
            data = data.reshape(np_, -1)
        if is_label and label_name is None:
            label_name = name
        point_data[name] = data
    regions = []
    table = root.find("LabelTable")
    if table is not None and label_name is not None:
        labels = np.asarray(point_data[label_name]).astype(np.int64)
        for label in table.findall("Label"):
            key = int(label.get("Key", label.get("Index", "0")))
            regions.append(
                Region(
                    (label.text or "").strip(),
                    "point",
                    np.flatnonzero(labels == key),
                    tag=key,
                )
            )
    return Mesh(
        points, cells, point_data=point_data, field_data=field_data, regions=regions
    )


def _array(intent, dtype_name, shape, raw, name, extra):
    dims = "".join(f' Dim{d}="{n}"' for d, n in enumerate(shape))
    head = (
        f'<DataArray Intent="{intent}" DataType="{dtype_name}" '
        f'ArrayIndexingOrder="RowMajorOrder" Dimensionality="{len(shape)}"{dims}'
        ' Encoding="Base64Binary" Endian="LittleEndian" ExternalFileName="" '
        'ExternalFileOffset="">\n'
    )
    meta = (
        "<MetaData/>\n"
        if not name
        else "<MetaData>\n<MD>\n<Name><![CDATA[Name]]></Name>\n<Value><![CDATA["
        + name
        + "]]></Value>\n</MD>\n</MetaData>\n"
    )
    data = base64.b64encode(raw).decode("ascii")
    return head + meta + extra + f"<Data>{data}</Data>\n</DataArray>\n"


def _text(mesh, key, default):
    value = mesh.field_data.get(key)
    if value is None:
        return default
    return np.asarray(value, dtype=np.uint8).tobytes().decode("utf-8", "replace")


def write(filename, mesh):
    if len(mesh.cells) > 1 or (
        len(mesh.cells) == 1 and mesh.cells[0].type != "triangle"
    ):
        raise WriteError("meshio++: gifti: a GIFTI surface is one triangle block")
    pts = np.asarray(mesh.points, dtype=np.float64)
    np_ = len(pts)
    full = np.zeros((np_, 3))
    full[:, : min(3, pts.shape[1])] = pts[:, :3]
    matrix = np.eye(4).reshape(-1)
    t = mesh.field_data.get("gifti:transform")
    if t is not None and np.asarray(t).size == 16:
        matrix = np.asarray(t, dtype=np.float64).reshape(-1)
    xform = (
        "<CoordinateSystemTransformMatrix>\n<DataSpace><![CDATA["
        + _text(mesh, "gifti:data_space", "NIFTI_XFORM_UNKNOWN")
        + "]]></DataSpace>\n<TransformedSpace><![CDATA["
        + _text(mesh, "gifti:transformed_space", "NIFTI_XFORM_UNKNOWN")
        + "]]></TransformedSpace>\n<MatrixData>"
        + " ".join("%.17g" % v for v in matrix)
        + "</MatrixData>\n</CoordinateSystemTransformMatrix>\n"
    )
    arrays = [
        _array(
            "NIFTI_INTENT_POINTSET",
            "NIFTI_TYPE_FLOAT32",
            (np_, 3),
            full.astype("<f4").tobytes(),
            "",
            xform,
        )
    ]
    if mesh.cells:
        tri = np.asarray(mesh.cells[0].data, dtype=np.int64)
        arrays.append(
            _array(
                "NIFTI_INTENT_TRIANGLE",
                "NIFTI_TYPE_INT32",
                (len(tri), 3),
                tri.astype("<i4").tobytes(),
                "",
                "",
            )
        )
    for name in sorted(mesh.point_data):
        a = np.asarray(mesh.point_data[name])
        if a.ndim == 0 or a.shape[0] != np_ or a.dtype.kind not in "fiu":
            continue
        width = a.size // np_ if np_ else 0
        floating = a.dtype.kind == "f"
        raw = a.reshape(-1).astype("<f4" if floating else "<i4").tobytes()
        if name == "gifti:label":
            intent = "NIFTI_INTENT_LABEL"
        elif width > 1 or not floating:
            intent = "NIFTI_INTENT_NONE"
        else:
            intent = "NIFTI_INTENT_SHAPE"
        arrays.append(
            _array(
                intent,
                "NIFTI_TYPE_FLOAT32" if floating else "NIFTI_TYPE_INT32",
                (np_, width) if width > 1 else (np_,),
                raw,
                name,
                "",
            )
        )
    labels = "".join(
        f'<Label Key="{int(r.tag)}" Red="0" Green="0" Blue="0" Alpha="0">'
        f"<![CDATA[{r.name}]]></Label>\n"
        for r in sorted(
            (r for r in mesh.regions if r.kind == "point" and r.tag >= 0),
            key=lambda r: (r.name, r.dim, r.tag),
        )
    )
    out = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<!DOCTYPE GIFTI SYSTEM "http://www.nitrc.org/frs/download.php/115/gifti.dtd">\n'
        f'<GIFTI Version="1.0" NumberOfDataArrays="{len(arrays)}">\n<MetaData/>\n'
        + (
            "<LabelTable/>\n"
            if not labels
            else "<LabelTable>\n" + labels + "</LabelTable>\n"
        )
        + "".join(arrays)
        + "</GIFTI>\n"
    )
    with open(filename, "wb") as f:
        f.write(out.encode("utf-8"))
