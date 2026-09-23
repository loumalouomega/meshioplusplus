#!/usr/bin/env python3
"""Regenerate the FEBio ``.feb`` fixtures under ``tests/python/meshes/febio/``.

The files are written here, not by meshio++, so the reader is never only checked
against its own writer. Each follows the syntax FEBio's own parsers accept for
its spec version (``FEBioXML/FEBioGeometrySection.cpp`` for 2.5,
``FEBioMeshSection.cpp`` for 3.0, ``FEBioMeshSection4.cpp`` for 4.0), and every
element lists its nodes in FEBio's numbering, from the shape functions in
``FECore/FESolidElementShape.cpp``:

* tet10, penta15, pyra13, hex20, tri6 and quad8: corners, then the mid-edge
  nodes in the same edge order as VTK;
* hex27: corners, the bottom, top and vertical mid-edges, then the mid-height
  face centres at ``s=-1, r=+1, s=+1, r=-1``, the bottom and top centres and
  the body centre (VTK puts the mid-height face centres at ``r=-1, r=+1, s=-1,
  s=+1``).

Three files:

* ``all_elements_v40.feb`` (spec 4.0): one element of every solid type, a
  shell, sparse node ids across two ``<Nodes>`` blocks, node sets with ranges,
  an element set, a surface on solid faces, a surface off them, an edge, a
  discrete set and node/element ``<MeshData>``. FEBio 4.12 reads it (checked
  with a ``<Control>`` section added).
* ``block_v30.feb`` (spec 3.0): two hex8 blocks with children-style sets.
* ``block_v25.feb`` (spec 2.5): the same in ``<Geometry>``, with ``mat=`` ids.

    python tools/gen_febio_fixtures.py
"""

import pathlib

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "febio"
)

HEX = [
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
]
HEX_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
HEX_EDGES += [(0, 4), (1, 5), (2, 6), (3, 7)]
TET = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
TET_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
WEDGE = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)]
WEDGE_EDGES = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]
PYRA = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0.5, 0.5, 1)]
PYRA_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)]


def mid(points):
    return tuple(sum(c) / len(points) for c in zip(*points))


class Nodes:
    """Coordinates -> ids, in order of first use, starting at `first` and
    stepping by `step` (sparse ids are legal in spec 4.0)."""

    def __init__(self, first=1, step=1):
        self.ids = {}
        self.first = first
        self.step = step

    def __call__(self, xyz):
        key = tuple(float(v) for v in xyz)
        if key not in self.ids:
            self.ids[key] = self.first + self.step * len(self.ids)
        return self.ids[key]


def shifted(corners, dx):
    return [(x + dx, y, z) for x, y, z in corners]


def quadratic(nodes, corners, edges, extra=()):
    ids = [nodes(p) for p in corners]
    ids += [nodes(mid([corners[a], corners[b]])) for a, b in edges]
    ids += [nodes(p) for p in extra]
    return ids


def hex27(nodes, corners):
    """FEBio hex27: corners, 12 edges, then face centres y-, x+, y+, x-, z-, z+,
    and the body centre."""
    c = corners
    faces = [
        (0, 1, 5, 4),
        (1, 2, 6, 5),
        (2, 3, 7, 6),
        (3, 0, 4, 7),
        (0, 1, 2, 3),
        (4, 5, 6, 7),
    ]
    extra = [mid([c[i] for i in f]) for f in faces] + [mid(c)]
    return quadratic(nodes, c, HEX_EDGES, extra)


def xyz_text(nodes):
    return {i: ",".join(f"{v:g}" for v in key) for key, i in nodes.ids.items()}


def all_elements_v40():
    solid = Nodes(first=1)
    blocks = [
        ("tet4", [[solid(p) for p in shifted(TET, 0)]]),
        ("tet10", [quadratic(solid, shifted(TET, 2), TET_EDGES)]),
        ("penta6", [[solid(p) for p in shifted(WEDGE, 4)]]),
        ("penta15", [quadratic(solid, shifted(WEDGE, 6), WEDGE_EDGES)]),
        ("pyra5", [[solid(p) for p in shifted(PYRA, 8)]]),
        ("pyra13", [quadratic(solid, shifted(PYRA, 10), PYRA_EDGES)]),
        ("hex8", [[solid(p) for p in shifted(HEX, 12)]]),
        ("hex20", [quadratic(solid, shifted(HEX, 14), HEX_EDGES)]),
        ("hex27", [hex27(solid, shifted(HEX, 16))]),
    ]
    # A shell plate in a second <Nodes> block with sparse ids.
    shell = Nodes(first=1001, step=2)
    plate = [(0, -2, 0), (1, -2, 0), (1, -1, 0), (0, -1, 0)]
    blocks.append(("quad4", [[shell(p) for p in plate]]))
    lines = [
        '<?xml version="1.0" encoding="ISO-8859-1"?>',
        '<febio_spec version="4.0">',
    ]
    lines.append('\t<Module type="solid"/>')
    lines.append("\t<Material>")
    for k, (etype, _) in enumerate(blocks, 1):
        lines.append(
            f'\t\t<material id="{k}" name="mat_{etype}" type="isotropic elastic">'
        )
        lines.append("\t\t\t<E>1</E>\n\t\t\t<v>0.3</v>\n\t\t</material>")
    lines.append("\t</Material>")
    lines.append("\t<Mesh>")
    for name, nodes in (("SolidNodes", solid), ("ShellNodes", shell)):
        lines.append(f'\t\t<Nodes name="{name}">')
        lines += [
            f'\t\t\t<node id="{i}">{t}</node>' for i, t in xyz_text(nodes).items()
        ]
        lines.append("\t\t</Nodes>")
    elem_id = 0
    for etype, rows in blocks:
        lines.append(f'\t\t<Elements type="{etype}" name="part_{etype}">')
        for row in rows:
            elem_id += 1
            lines.append(f'\t\t\t<elem id="{elem_id}">{",".join(map(str, row))}</elem>')
        lines.append("\t\t</Elements>")
    # hex8 bottom face (outward) and a triangle off every solid face.
    hex8 = blocks[6][1][0]
    lines.append(
        '\t\t<NodeSet name="hex_bottom">' + ",".join(map(str, hex8[:4])) + "</NodeSet>"
    )
    lines.append('\t\t<NodeSet name="first_ten">1:10</NodeSet>')
    lines.append('\t\t<ElementSet name="quadratic">2,4,6,8,9</ElementSet>')
    lines.append('\t\t<Surface name="hex_base">')
    lines.append(f'\t\t\t<quad4 id="1">{hex8[0]},{hex8[3]},{hex8[2]},{hex8[1]}</quad4>')
    lines.append("\t\t</Surface>")
    lines.append('\t\t<Surface name="floating">')
    lines.append(f'\t\t\t<tri3 id="1">{hex8[0]},{hex8[1]},{hex8[6]}</tri3>')
    lines.append("\t\t</Surface>")
    lines.append('\t\t<Edge name="hex_edge">')
    lines.append(f'\t\t\t<line2 id="1">{hex8[0]},{hex8[1]}</line2>')
    lines.append("\t\t</Edge>")
    lines.append('\t\t<DiscreteSet name="springs">')
    lines.append(f"\t\t\t<delem>{hex8[6]},{blocks[8][1][0][0]}</delem>")
    lines.append("\t\t</DiscreteSet>")
    lines.append("\t</Mesh>")
    lines.append("\t<MeshDomains>")
    for etype, _ in blocks:
        if etype == "quad4":
            lines.append(f'\t\t<ShellDomain name="part_{etype}" mat="mat_{etype}">')
            lines.append(
                "\t\t\t<shell_thickness>0.1</shell_thickness>\n\t\t</ShellDomain>"
            )
        else:
            lines.append(f'\t\t<SolidDomain name="part_{etype}" mat="mat_{etype}"/>')
    lines.append("\t</MeshDomains>")
    lines.append("\t<MeshData>")
    lines.append(
        '\t\t<NodeData name="temperature" node_set="hex_bottom" data_type="scalar">'
    )
    lines += [f'\t\t\t<node lid="{k}">{10 * k}</node>' for k in range(1, 5)]
    lines.append("\t\t</NodeData>")
    # FEBio wants element data on one domain's elements: the hex8 block's set.
    lines.append('\t\t<ElementData name="fiber" elem_set="part_hex8" data_type="vec3">')
    lines.append('\t\t\t<e lid="1">1,0,0</e>')
    lines.append("\t\t</ElementData>")
    lines.append("\t</MeshData>")
    lines.append("</febio_spec>")
    return "\n".join(lines) + "\n"


def two_blocks(version):
    """Two hex8 blocks side by side, with sets in the children form."""
    nodes = Nodes()
    left = [nodes(p) for p in HEX]
    right = [nodes(p) for p in shifted(HEX, 1)]
    spec25 = version == "2.5"
    section = "Geometry" if spec25 else "Mesh"
    lines = [
        '<?xml version="1.0" encoding="ISO-8859-1"?>',
        f'<febio_spec version="{version}">',
    ]
    lines.append('\t<Module type="solid"/>')
    lines.append("\t<Material>")
    lines.append(
        '\t\t<material id="1" name="soft" type="isotropic elastic"><E>1</E><v>0.3</v></material>'
    )
    lines.append(
        '\t\t<material id="2" name="stiff" type="isotropic elastic"><E>9</E><v>0.3</v></material>'
    )
    lines.append("\t</Material>")
    lines.append(f"\t<{section}>")
    lines.append('\t\t<Nodes name="all">')
    lines += [f'\t\t\t<node id="{i}">{t}</node>' for i, t in xyz_text(nodes).items()]
    lines.append("\t\t</Nodes>")
    for eid, name, mat, row in ((1, "left", 1, left), (2, "right", 2, right)):
        mat_attr = f' mat="{mat}"' if spec25 else ""
        lines.append(f'\t\t<Elements type="hex8" name="{name}"{mat_attr}>')
        lines.append(f'\t\t\t<elem id="{eid}">{",".join(map(str, row))}</elem>')
        lines.append("\t\t</Elements>")
    child = "node"
    lines.append('\t\t<NodeSet name="fixed">')
    lines += [f'\t\t\t<{child} id="{n}"/>' for n in left[:4]]
    lines.append("\t\t</NodeSet>")
    lines.append('\t\t<ElementSet name="both">')
    lines += [f'\t\t\t<elem id="{e}"/>' for e in (1, 2)]
    lines.append("\t\t</ElementSet>")
    lines.append('\t\t<Surface name="top">')
    for k, row in enumerate((left, right), 1):
        lines.append(
            f'\t\t\t<quad4 id="{k}">{row[4]},{row[5]},{row[6]},{row[7]}</quad4>'
        )
    lines.append("\t\t</Surface>")
    lines.append(f"\t</{section}>")
    if not spec25:
        lines.append("\t<MeshDomains>")
        lines.append('\t\t<SolidDomain name="left" mat="soft"/>')
        lines.append('\t\t<SolidDomain name="right" mat="stiff"/>')
        lines.append("\t</MeshDomains>")
    lines.append("\t<MeshData>")
    tag = "elem" if spec25 else "e"
    lines.append('\t\t<ElementData name="thickness" elem_set="both" datatype="scalar">')
    lines += [f'\t\t\t<{tag} lid="{k}">{0.5 * k}</{tag}>' for k in (1, 2)]
    lines.append("\t\t</ElementData>")
    lines.append("\t</MeshData>")
    lines.append("</febio_spec>")
    return "\n".join(lines) + "\n"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, text in (
        ("all_elements_v40.feb", all_elements_v40()),
        ("block_v30.feb", two_blocks("3.0")),
        ("block_v25.feb", two_blocks("2.5")),
    ):
        with open(OUT / name, "w", newline="\n") as f:
            f.write(text)
        print(f"wrote {OUT / name}")


if __name__ == "__main__":
    main()
