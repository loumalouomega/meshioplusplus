"""Named groups of mesh entities — the Python side of the core ``Region`` model.

A *region* is a named group of points, cells or cell **facets**: what gmsh calls
a physical group, Exodus a block / node set / side set, Abaqus an ``*NSET`` /
``*ELSET`` / ``*SURFACE``, MED a family or group, UNV a group, Ansys a component,
OpenFOAM a boundary patch and Kratos a SubModelPart. Regions are first-class
members of the C++ core (``src/cpp/include/meshioplusplus/region.hpp``), so they
cross the pybind11 boundary natively and are visible from the C API, Fortran,
WebAssembly and the native CLI.

Three kinds, spelled by the ``kind`` string:

``"point"``
    ``entries`` are point indices, shape ``(n,)``.
``"cell"``
    ``entries`` are **global** cell indices, shape ``(n,)``, numbered
    *block-major*: block 0's cells first, then block 1's, and so on (the same
    numbering :func:`meshioplusplus.partition_labels` uses).
``"side"``
    ``entries`` are ``(global cell index, local facet index)`` pairs, shape
    ``(n, 2)``. Facets are numbered as the core's ``detail/cell_faces.hpp``
    numbers the faces of a 3-D cell and ``detail/cell_edges.hpp`` the edges of a
    2-D one. This kind has no equivalent in the legacy ``point_sets`` /
    ``cell_sets`` model.

Legacy compatibility
--------------------
``Mesh.point_sets`` and ``Mesh.cell_sets`` keep working exactly as they always
have. They are now **views** over the mesh's regions: reading one materializes
the historical shape (a flat index array for ``point_sets``, a per-block list of
local index arrays for ``cell_sets``), and writing one creates or replaces the
matching regions. No deprecation, no behaviour change.

Two deliberate accommodations make that true rather than nearly true:

* **The escape hatch.** Formats stash things in ``cell_sets`` that are not cell
  indices at all — gmsh's ``gmsh:bounding_entities`` holds entity tags, negative
  values included. Such an entry cannot be a region (regions are sorted,
  de-duplicated, non-negative indices), so it is kept verbatim in a private
  passthrough dict and returned unchanged. :func:`cell_set_is_index_data` is the
  predicate, shared with :mod:`meshioplusplus._split`.
* **``None`` blocks become empty arrays.** A ``cell_sets`` value may carry
  ``None`` for a block with no members (gmsh's reader builds ``[None] * n``).
  Global indices cannot express "absent" separately from "empty", so a
  round-trip through a region normalizes ``None`` to an empty array. Every
  consumer in the tree already treats the two identically — and the Abaqus
  writer's ``len(v[ic])`` would in fact raise on a genuine ``None``.
"""

from __future__ import annotations

from typing import Union

import numpy as np

#: The region kinds, in the core's enumerator order (``RegionKind``).
KINDS = ("point", "cell", "side")


class Region:
    """A named group of points, cells or cell facets.

    Parameters
    ----------
    name :
        the group's name, as the originating format spelled it.
    kind :
        one of ``"point"``, ``"cell"``, ``"side"`` (see the module docstring).
    entries :
        int64 indices — shape ``(n,)`` for point/cell, ``(n, 2)`` for side.
        Sorted and de-duplicated on construction, so two regions describing the
        same group compare equal.
    dim :
        topological dimension the group was declared for, or ``-1`` when the
        format does not say (gmsh physical groups are per-dimension).
    tag :
        the format-native integer id (gmsh physical tag, MED family id, Exodus
        set id), or ``-1`` when there is none.
    """

    __slots__ = ("name", "kind", "dim", "tag", "entries")

    def __init__(self, name, kind, entries, dim: int = -1, tag: int = -1):
        if kind not in KINDS:
            raise ValueError(
                f"region '{name}': unknown kind '{kind}' "
                f"(expected one of {', '.join(KINDS)})"
            )
        self.name = str(name)
        self.kind = kind
        self.dim = int(dim)
        self.tag = int(tag)
        self.entries = _canonical_entries(entries, kind)

    def __len__(self):
        return len(self.entries)

    def __repr__(self):
        return (
            f"<meshio++ Region {self.name!r} kind={self.kind} "
            f"dim={self.dim} tag={self.tag} entries={len(self.entries)}>"
        )

    def __eq__(self, other):
        if not isinstance(other, Region):
            return NotImplemented
        return (
            self.name == other.name
            and self.kind == other.kind
            and self.dim == other.dim
            and self.tag == other.tag
            and self.entries.shape == other.entries.shape
            and bool(np.array_equal(self.entries, other.entries))
        )

    @property
    def key(self):
        """``(kind, name, dim, tag)`` — the identity two regions replace on."""
        return (KINDS.index(self.kind), self.name, self.dim, self.tag)

    def copy(self) -> Region:
        """A deep copy (the entries array is copied, not shared)."""
        return Region(self.name, self.kind, self.entries.copy(), self.dim, self.tag)


def _canonical_entries(entries, kind):
    """Sort, de-duplicate and shape an entry array for ``kind``."""
    arr = np.asarray(entries)
    if arr.size == 0:
        return np.empty((0, 2) if kind == "side" else (0,), dtype=np.int64)
    arr = arr.astype(np.int64, copy=False)
    if kind == "side":
        arr = arr.reshape(-1, 2)
        # Lexicographic on (cell, facet); np.unique on a structured view is the
        # numpy idiom for "unique rows", and it sorts as a side effect.
        return np.unique(arr, axis=0)
    return np.unique(arr.reshape(-1))


# --------------------------------------------------------------------------- #
# global (block-major) cell index                                             #
# --------------------------------------------------------------------------- #
def block_bases(cells) -> np.ndarray:
    """Prefix sums of the per-block cell counts.

    Entry ``b`` is the global index of block ``b``'s first cell; the last entry
    is the total cell count. Every block contributes, including empty and ragged
    ones, so the numbering never skips — this is the Python twin of the core's
    ``detail/cell_index.hpp``.
    """
    bases = np.zeros(len(cells) + 1, dtype=np.int64)
    for b, cb in enumerate(cells):
        bases[b + 1] = bases[b] + len(cb.data)
    return bases


def cell_set_is_index_data(cells, blocks) -> bool:
    """Whether a ``cell_sets`` value is per-block cell-index lists.

    Formats stash other per-block side data in ``cell_sets`` — gmsh's
    ``gmsh:bounding_entities`` holds entity tags, not cell indices. Such an entry
    cannot become a region; it takes the passthrough escape hatch instead.

    The test is **non-negative integers**, not "integers in range". A negative
    value is the reliable tell: entity tags carry an orientation sign, and in the
    bundled ``example.msh`` 475 of the 476 populated blocks contain one, so the
    whole entry is rejected on that alone. Requiring every index to be *in range*
    for its block would instead reject a genuine — if sloppily built — cell set
    whose indices overshoot, and meshio has always accepted those and quietly
    dropped the overshooting entries on a round-trip. :func:`blocks_to_global`
    reproduces that by skipping them.
    """
    if isinstance(blocks, np.ndarray) or not hasattr(blocks, "__len__"):
        return False
    for b in range(len(cells)):
        if b >= len(blocks) or blocks[b] is None:
            continue
        try:
            idx = np.asarray(blocks[b]).reshape(-1)
        except (ValueError, TypeError):
            return False
        if len(idx) == 0:
            continue
        if idx.dtype.kind not in "iu":
            return False
        if idx.min() < 0:
            return False
    return True


def blocks_to_global(cells, blocks) -> np.ndarray:
    """Per-block local cell indices -> global (block-major) indices.

    Indices past the end of their block are dropped rather than producing a
    global index belonging to the *next* block — see
    :func:`cell_set_is_index_data` for why such entries are accepted at all.
    """
    bases = block_bases(cells)
    out = []
    for b in range(len(cells)):
        if b >= len(blocks) or blocks[b] is None:
            continue
        idx = np.asarray(blocks[b], dtype=np.int64).reshape(-1)
        idx = idx[idx < len(cells[b].data)]
        if len(idx):
            out.append(idx + bases[b])
    if not out:
        return np.empty(0, dtype=np.int64)
    return np.concatenate(out)


def global_to_blocks(cells, entries) -> list:
    """Global (block-major) cell indices -> one local index array per block.

    Always returns exactly ``len(cells)`` arrays, so the result has the shape
    every existing ``cell_sets`` consumer expects.
    """
    bases = block_bases(cells)
    idx = np.asarray(entries, dtype=np.int64).reshape(-1)
    out = []
    for b in range(len(cells)):
        lo, hi = bases[b], bases[b + 1]
        sel = idx[(idx >= lo) & (idx < hi)] - lo
        out.append(sel.astype(np.int64))
    return out


# --------------------------------------------------------------------------- #
# the compat views                                                            #
# --------------------------------------------------------------------------- #
class _SetsView(dict):
    """Base of the ``point_sets`` / ``cell_sets`` views over ``mesh.regions``.

    Deliberately a ``dict`` **subclass**, not a bare ``MutableMapping``: plenty
    of code — including ``numpy.testing.assert_equal``, which this repo's own
    tests use on ``mesh.cell_sets`` — branches on ``isinstance(x, dict)``, so
    anything less would be an observable break of the legacy API.

    Each access to ``mesh.point_sets`` / ``mesh.cell_sets`` builds a fresh view
    holding a *snapshot* of the current regions; the mutating methods are
    overridden to write straight through to ``mesh.regions``, so
    ``mesh.cell_sets[name] = value`` and ``mesh.cell_sets.update(...)`` behave
    exactly as they always did. (Note that ``dict.update`` does not route
    through ``__setitem__``, which is why every mutator is overridden here
    rather than inherited.)
    """

    _kind = "point"

    def __init__(self, mesh):
        super().__init__(self._materialize(mesh))
        self._mesh = mesh

    @staticmethod
    def _materialize(mesh):  # pragma: no cover - overridden in both subclasses
        return {}

    # -- the region slice this view exposes --------------------------------- #
    def _find_index(self, name):
        """Index into ``mesh.regions`` of the first region of this kind and name.

        Index-based rather than object-based: ``list.index`` compares by value,
        and a mesh may legitimately hold two regions with the same name and kind
        but different ``dim``/``tag`` (gmsh physical groups of two dimensions).
        """
        for i, r in enumerate(self._mesh.regions):
            if r.kind == self._kind and r.name == name:
                return i
        return -1

    def _find(self, name):
        i = self._find_index(name)
        return self._mesh.regions[i] if i >= 0 else None

    def _store(self, name, region):
        """Append the region, or replace the first one with this name and kind."""
        i = self._find_index(name)
        if i < 0:
            self._mesh.regions.append(region)
        else:
            self._mesh.regions[i] = region

    def _drop(self, name):
        self._mesh.regions[:] = [
            r
            for r in self._mesh.regions
            if not (r.kind == self._kind and r.name == name)
        ]

    # -- mutators: write through, then update the snapshot ------------------- #
    def __delitem__(self, name):
        if name not in self:
            raise KeyError(name)
        self._drop(name)
        super().__delitem__(name)

    def update(self, other=(), **kwargs):
        items = other.items() if hasattr(other, "items") else other
        for name, value in items:
            self[name] = value
        for name, value in kwargs.items():
            self[name] = value

    def setdefault(self, name, default=None):
        if name not in self:
            self[name] = default
        return self[name]

    def pop(self, name, *default):
        if name not in self:
            if default:
                return default[0]
            raise KeyError(name)
        value = self[name]
        del self[name]
        return value

    def popitem(self):
        if not self:
            raise KeyError("popitem(): dictionary is empty")
        name = next(reversed(self))
        value = self[name]
        del self[name]
        return name, value

    def clear(self):
        for name in list(self):
            del self[name]


class PointSetsView(_SetsView):
    """``mesh.point_sets``: a dict of ``name -> int64 point-index array``."""

    _kind = "point"

    @staticmethod
    def _materialize(mesh):
        return {r.name: r.entries for r in mesh.regions if r.kind == "point"}

    def __setitem__(self, name, value):
        existing = self._find(name)
        # Setting through the view must not silently discard a dim/tag the
        # originating format supplied.
        dim = existing.dim if existing is not None else -1
        tag = existing.tag if existing is not None else -1
        region = Region(name, "point", value, dim=dim, tag=tag)
        self._store(name, region)
        super().__setitem__(name, region.entries)


class CellSetsView(_SetsView):
    """``mesh.cell_sets``: a dict of ``name -> list of per-block index arrays``.

    The list always has exactly ``len(mesh.cells)`` entries, materialized from
    the region's global indices against the mesh's *current* blocks. Entries
    that are not cell-index data (see :func:`cell_set_is_index_data`) live in the
    mesh's passthrough dict and are returned verbatim.
    """

    _kind = "cell"

    @staticmethod
    def _materialize(mesh):
        out = {
            r.name: global_to_blocks(mesh.cells, r.entries)
            for r in mesh.regions
            if r.kind == "cell"
        }
        out.update(mesh._extra_cell_sets)
        return out

    def __setitem__(self, name, value):
        extra = self._mesh._extra_cell_sets
        if not cell_set_is_index_data(self._mesh.cells, value):
            # Not cell indices (gmsh's entity tags, a bare array, ...) — keep it
            # verbatim so `mesh.cell_sets` stays lossless.
            self._drop(name)
            extra[name] = value
            super().__setitem__(name, value)
            return
        extra.pop(name, None)
        existing = self._find(name)
        dim = existing.dim if existing is not None else -1
        tag = existing.tag if existing is not None else -1
        region = Region(
            name,
            "cell",
            blocks_to_global(self._mesh.cells, value),
            dim=dim,
            tag=tag,
        )
        self._store(name, region)
        super().__setitem__(name, global_to_blocks(self._mesh.cells, region.entries))

    def __delitem__(self, name):
        extra = self._mesh._extra_cell_sets
        if name in extra:
            del extra[name]
            dict.__delitem__(self, name)
            return
        super().__delitem__(name)


# --------------------------------------------------------------------------- #
# helpers used by the format shims and operations                             #
# --------------------------------------------------------------------------- #
def regions_of_kind(mesh, kind) -> list:
    """Every region of ``mesh`` with the given kind, in mesh order."""
    return [r for r in getattr(mesh, "regions", []) if r.kind == kind]


def replace_regions(mesh, regions) -> None:
    """Replace ``mesh``'s whole region list (used by the operation shims)."""
    mesh.regions[:] = list(regions)


def side_entries(pairs) -> np.ndarray:
    """Build a canonical ``(n, 2)`` side-entry array from ``(cell, facet)`` pairs."""
    return _canonical_entries(np.asarray(list(pairs), dtype=np.int64), "side")


def format_regions(mesh) -> Union[str, None]:
    """One human-readable line per kind, for ``Mesh.__repr__`` / the CLI."""
    regions = getattr(mesh, "regions", None)
    if not regions:
        return None
    lines = []
    for kind, label in (
        ("point", "Point sets"),
        ("cell", "Cell sets"),
        ("side", "Side sets"),
    ):
        names = [r.name for r in regions if r.kind == kind]
        if names:
            lines.append(f"  {label}: {', '.join(names)}")
    return "\n".join(lines) if lines else None
