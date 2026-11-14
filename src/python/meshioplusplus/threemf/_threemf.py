"""3MF packages: the Python reference of ``formats/threemf.cpp`` (roadmap 2.4.5). See
that header for what is read and written; both engines write the same bytes."""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET
import zipfile

import numpy as np

from .._exceptions import ReadError, WriteError
from .._mesh import Mesh
from .._regions import Region
from .._units import LENGTH_UNIT_KEY, length_unit
from .._zip import ZipWriter

COLOR = "3mf:color"
MODEL_TYPE = "http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"
UNITS = [
    ("micron", 1e-6),
    ("millimeter", 1e-3),
    ("centimeter", 1e-2),
    ("inch", 0.0254),
    ("foot", 0.3048),
    ("meter", 1.0),
]
_IDENTITY = (1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0)


def _local(tag):
    return tag.rsplit("}", 1)[-1].rsplit(":", 1)[-1]


def _child(node, name):
    for c in node:
        if _local(c.tag) == name:
            return c
    return None


def _children(node, name):
    return [c for c in (node if node is not None else []) if _local(c.tag) == name]


def _number(node, attr):
    try:
        v = float(node.get(attr, ""))
    except ValueError:
        v = math.nan
    if not math.isfinite(v):
        raise ReadError(f"meshio++: 3mf: <{_local(node.tag)}> has no number '{attr}'")
    return v


def _transform(text):
    if not text:
        return _IDENTITY
    try:
        m = tuple(float(v) for v in text.split()[:12])
    except ValueError:
        m = ()
    if len(m) != 12:
        raise ReadError("meshio++: 3mf: a transform needs twelve numbers")
    return m


def _compose(a, b):
    out = []
    for r in range(4):
        for c in range(3):
            v = b[9 + c] if r == 3 else 0.0
            for k in range(3):
                v += a[r * 3 + k] * b[k * 3 + c]
            out.append(v)
    return tuple(out)


def _colour(text):
    rgba = [0, 0, 0, 255]
    if len(text) in (7, 9):
        for k in range(0, len(text) - 1, 2):
            rgba[k // 2] = int(text[1 + k : 3 + k], 16)
    return rgba


def _int(node, attr, default):
    try:
        return int(node.get(attr))
    except (TypeError, ValueError):
        return default


def _instance(objects, materials, oid, m, out, depth):
    if depth > 64:
        raise ReadError("meshio++: 3mf: components nest deeper than 64 (a cycle?)")
    if oid not in objects:
        raise ReadError(
            f"meshio++: 3mf: the build names object {oid}, which the model does not have"
        )
    obj = objects[oid]
    mesh = _child(obj, "mesh")
    if mesh is not None:
        base = len(out["points"])
        verts = _children(_child(mesh, "vertices"), "vertex")
        for v in verts:
            x, y, z = _number(v, "x"), _number(v, "y"), _number(v, "z")
            out["points"].append(
                (
                    x * m[0] + y * m[3] + z * m[6] + m[9],
                    x * m[1] + y * m[4] + z * m[7] + m[10],
                    x * m[2] + y * m[5] + z * m[8] + m[11],
                )
            )
        obj_pid = _int(obj, "pid", -1)
        obj_index = _int(obj, "pindex", 0)
        for t in _children(_child(mesh, "triangles"), "triangle"):
            row = []
            for a in ("v1", "v2", "v3"):
                v = _int(t, a, -1)
                if not 0 <= v < len(verts):
                    raise ReadError(
                        f"meshio++: 3mf: a triangle names vertex {v} of {len(verts)}"
                    )
                row.append(base + v)
            out["triangles"].append(row)
            pid = _int(t, "pid", obj_pid)
            index = _int(t, "p1", 0 if t.get("pid") is not None else obj_index)
            colours = materials.get(pid)
            if colours is not None and 0 <= index < len(colours):
                out["colours"].append(colours[index])
                out["any"] = True
            else:
                out["colours"].append([0, 0, 0, 0])
        return
    comps = _child(obj, "components")
    for c in _children(comps, "component"):
        _instance(
            objects,
            materials,
            _int(c, "objectid", -1),
            _compose(_transform(c.get("transform")), m),
            out,
            depth + 1,
        )


def read(filename):
    try:
        archive = zipfile.ZipFile(str(filename))
    except (zipfile.BadZipFile, OSError) as exc:
        raise ReadError(
            f"meshio++: 3mf: '{filename}' is not a ZIP package: {exc}"
        ) from None
    with archive:
        names = set(archive.namelist())
        model_name = "3D/3dmodel.model"
        if "_rels/.rels" in names:
            try:
                rels = ET.fromstring(archive.read("_rels/.rels"))
                for r in rels:
                    if r.get("Type") == MODEL_TYPE:
                        model_name = r.get("Target", "").lstrip("/")
                        break
            except ET.ParseError:
                pass
        if model_name not in names:
            raise ReadError(
                f"meshio++: 3mf: the package has no model part '{model_name}'"
            )
        try:
            model = ET.fromstring(archive.read(model_name))
        except ET.ParseError:
            raise ReadError("meshio++: 3mf: the model part is not XML") from None
    if _local(model.tag) != "model":
        raise ReadError("meshio++: 3mf: the model part has no <model> root")
    required = model.get("requiredextensions", "")
    if required.strip():
        raise ReadError(
            f"meshio++: 3mf: the model requires the extensions '{required}', which "
            "meshio++ does not read (the core specification only)"
        )
    objects, materials = {}, {}
    resources = _child(model, "resources")
    for r in resources if resources is not None else []:
        name = _local(r.tag)
        if name == "object":
            objects[_int(r, "id", -1)] = r
        elif name == "basematerials":
            materials[_int(r, "id", -1)] = [
                _colour(b.get("displaycolor", "")) for b in _children(r, "base")
            ]
    out = {"points": [], "triangles": [], "colours": [], "any": False}
    items, seen = [], {}
    for item in _children(_child(model, "build"), "item"):
        oid = _int(item, "objectid", -1)
        first = len(out["triangles"])
        _instance(objects, materials, oid, _transform(item.get("transform")), out, 0)
        name = objects[oid].get("name", "") if oid in objects else ""
        name = name or f"object_{oid}"
        n = seen.get(name, 0)
        seen[name] = n + 1
        if n:
            name = f"{name}_{n}"
        items.append((name, first, len(out["triangles"])))
    points = np.asarray(out["points"], dtype=np.float64).reshape(-1, 3)
    tris = np.asarray(out["triangles"], dtype=np.int64).reshape(-1, 3)
    cells = [("triangle", tris)] if len(tris) else []
    cell_data = {}
    if len(tris) and out["any"]:
        cell_data[COLOR] = [np.asarray(out["colours"], dtype=np.uint8)]
    regions = [Region(n, "cell", np.arange(a, b)) for n, a, b in items]
    unit = model.get("unit", "")
    metres = dict(UNITS).get(unit, 1e-3)
    return Mesh(
        points,
        cells,
        cell_data=cell_data,
        field_data={LENGTH_UNIT_KEY: np.array([metres])},
        regions=regions,
    )


def _escape(text):
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write(filename, mesh):
    tris = []
    for block in mesh.cells:
        if block.type != "triangle":
            raise WriteError(
                f"meshio++: 3mf: a '{block.type}' block; 3MF meshes are triangles"
            )
        tris.extend(np.asarray(block.data, dtype=np.int64).tolist())
    owner = [-1] * len(tris)
    names = []
    for r in sorted(
        (r for r in mesh.regions if r.kind == "cell"),
        key=lambda r: (r.name, r.dim, r.tag),
    ):
        used = False
        for c in np.asarray(r.entries, dtype=np.int64).reshape(-1):
            if 0 <= c < len(tris) and owner[c] < 0:
                owner[c] = len(names)
                used = True
        if used:
            names.append(r.name)
    rest = len(names)
    if -1 in owner or not tris:
        names.append("mesh")
    owner = [rest if o < 0 else o for o in owner]
    colours, colour_of = [], [-1] * len(tris)
    data = mesh.cell_data.get(COLOR)
    if data is not None and len(data) == len(mesh.cells):
        t = 0
        for block, arr in zip(mesh.cells, data):
            arr = np.asarray(arr).reshape(-1)
            if arr.size != len(block.data) * 4:
                raise WriteError(
                    "meshio++: 3mf: 3mf:color needs four uint8 values per triangle"
                )
            for row in arr.reshape(-1, 4).astype(np.int64) & 0xFF:
                c = tuple(int(v) for v in row)
                if c not in colours:
                    colours.append(c)
                colour_of[t] = colours.index(c)
                t += 1
    unit = "millimeter"
    metres = length_unit(mesh)
    for name, scale in UNITS:
        if metres is not None and abs(metres - scale) <= 1e-12 * scale:
            unit = name
    x = [
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<model unit="{unit}" xml:lang="en-US" '
        'xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">\n'
        '<metadata name="Application">meshio++</metadata>\n<resources>\n'
    ]
    next_id = 1
    material_id = -1
    if colours:
        material_id = next_id
        next_id += 1
        x.append(f'<basematerials id="{material_id}">\n')
        for k, c in enumerate(colours):
            x.append(
                f'<base name="color_{k}" displaycolor="#{c[0]:02X}{c[1]:02X}{c[2]:02X}{c[3]:02X}"/>\n'
            )
        x.append("</basematerials>\n")
    pts = np.asarray(mesh.points, dtype=np.float64)
    pd = pts.shape[1]
    ids = []
    for o, name in enumerate(names):
        oid = next_id
        next_id += 1
        ids.append(oid)
        x.append(
            f'<object id="{oid}" type="model" name="{_escape(name)}">\n<mesh>\n<vertices>\n'
        )
        local, order = {}, []
        for t, tri in enumerate(tris):
            if owner[t] == o:
                for v in tri:
                    if v not in local:
                        local[v] = len(order)
                        order.append(v)
        for v in order:
            coords = [pts[v, c] if c < pd else 0.0 for c in range(3)]
            x.append(
                "<vertex"
                + "".join(f' {a}="{"%.17g" % c}"' for a, c in zip("xyz", coords))
                + "/>\n"
            )
        x.append("</vertices>\n<triangles>\n")
        for t, tri in enumerate(tris):
            if owner[t] != o:
                continue
            row = f'<triangle v1="{local[tri[0]]}" v2="{local[tri[1]]}" v3="{local[tri[2]]}"'
            if material_id > 0:
                row += f' pid="{material_id}" p1="{colour_of[t]}"'
            x.append(row + "/>\n")
        x.append("</triangles>\n</mesh>\n</object>\n")
    x.append("</resources>\n<build>\n")
    x.extend(f'<item objectid="{i}"/>\n' for i in ids)
    x.append("</build>\n</model>\n")
    zip_ = ZipWriter()
    zip_.add(
        "[Content_Types].xml",
        '<?xml version="1.0" encoding="UTF-8"?>\n<Types '
        'xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" '
        'ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="model" '
        'ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>\n',
    )
    zip_.add(
        "_rels/.rels",
        '<?xml version="1.0" encoding="UTF-8"?>\n<Relationships '
        'xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Target="/3D/3dmodel.model" Id="rel0" Type="'
        + MODEL_TYPE
        + '"/></Relationships>\n',
    )
    zip_.add("3D/3dmodel.model", "".join(x))
    with open(filename, "wb") as f:
        f.write(zip_.finish())
