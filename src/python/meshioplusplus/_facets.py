"""Find the cell facet a file names by its nodes.

The Python twin of ``detail/facet_index.hpp``: LS-DYNA segments, FEBio
surfaces and Elmer boundary elements name a facet by its node list, and a
``side`` region names it by (global cell, local facet) -- faces for a 3-D cell,
edges for a 2-D one (``doc/regions.md``).
"""

from __future__ import annotations

from ._skin import _CELL_FACES
from ._surface import _CELL_EDGES


class FacetHit:
    """The owners of one facet: the first two in block-major order, and how many."""

    __slots__ = ("first", "second", "count")

    def __init__(self, owner):
        self.first = owner
        self.second = None
        self.count = 1


class FacetIndex:
    """Sorted-corner lookup of every facet of a mesh.

    ``solid_faces``: the faces of every 3-D cell. ``surface_edges``: the edges
    of every 2-D cell. ``surface_self``: a ``triangle``/``quad`` cell's own face
    as facet 0 (not a side-region facet; LS-DYNA shells and FEBio shell
    surfaces use it).
    """

    def __init__(self, mesh, solid_faces=True, surface_edges=True, surface_self=False):
        self._map = {}
        base = 0
        for block in mesh.cells:
            data = block.data
            faces = _CELL_FACES.get(block.type) if solid_faces else None
            edges = _CELL_EDGES.get(block.type) if surface_edges else None
            is_self = surface_self and block.type in ("triangle", "quad")
            if isinstance(data, list) or not (faces or edges or is_self):
                base += len(data)
                continue
            for r, row in enumerate(data.tolist()):
                g = base + r
                if faces:
                    for k, (_, ncorner, local) in enumerate(faces):
                        self._add(
                            tuple(sorted(row[i] for i in local[:ncorner])), (g, k)
                        )
                    continue
                if is_self:
                    self._add(tuple(sorted(row)), (g, 0))
                if edges:
                    for k, (_, _, local) in enumerate(edges):
                        self._add(tuple(sorted((row[local[0]], row[local[1]]))), (g, k))
            base += len(data)

    def _add(self, key, owner):
        hit = self._map.get(key)
        if hit is None:
            self._map[key] = FacetHit(owner)
            return
        if hit.count == 1:
            hit.second = owner
        hit.count += 1

    def find(self, corners):
        """The :class:`FacetHit` of the facet with these corner nodes, or ``None``."""
        return self._map.get(tuple(sorted(int(c) for c in corners)))

    def __len__(self):
        return len(self._map)


def facet_nodes(mesh, cell, facet):
    """``(facet_type, nodes)`` of one side facet, or ``None`` when it does not exist."""
    base = 0
    for block in mesh.cells:
        n = len(block.data)
        if cell < base + n:
            if isinstance(block.data, list) or facet < 0:
                return None
            row = block.data[cell - base]
            table = _CELL_FACES.get(block.type) or _CELL_EDGES.get(block.type)
            if not table or facet >= len(table):
                return None
            ftype, _, local = table[facet]
            return ftype, [int(row[i]) for i in local]
        base += n
    return None
