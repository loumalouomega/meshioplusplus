"""Set algebra and bookkeeping on named regions.

``edit_regions(mesh, edits)`` applies a list of edits, in order, to a copy of
the mesh's regions -- union, intersection, difference, rename, retag and
delete -- and leaves points, cells and data untouched. The semantics are
``operations/region_ops.hpp``'s, which the C++ core, the CLIs, the pipeline and
the flat bindings use.

This is implemented in Python rather than by a round trip through the core:
the edit touches only the region list, and handing the mesh to the core would
copy every array for nothing. ``tests/python/test_region_ops.py`` pins it to
``_core.edit_regions``.

An edit is a dict::

    {"op": "union", "inputs": ["inlet_a", "inlet_b"], "output": "inlet"}
    {"op": "difference", "inputs": ["skin", "inlet", "outlet"], "output": "wall"}
    {"op": "rename", "inputs": ["Surface 3"], "output": "Outlet"}
    {"op": "retag", "inputs": [{"name": "wall", "kind": "cell"}], "tag": 10}
    {"op": "delete", "inputs": ["scratch"]}

An input is a region name, or a dict ``{"name", "kind", "dim", "tag"}`` pinning
the rest of its identity; it must match exactly one region.
"""

from __future__ import annotations

import numpy as np

from ._regions import KINDS, Region

_PREFIX = "meshio++: edit_regions: "
_OPS = ("union", "intersection", "difference", "rename", "retag", "delete")


def _op_name(op):
    op = str(op)
    if op == "intersect":
        op = "intersection"
    if op not in _OPS:
        raise ValueError(
            f"{_PREFIX}unknown operation '{op}' (expected union, intersection, "
            "difference, rename, retag or delete)"
        )
    return op


def _selector(spec):
    if isinstance(spec, str):
        return {"name": spec, "kind": None, "dim": None, "tag": None}
    if isinstance(spec, Region):
        return {"name": spec.name, "kind": spec.kind, "dim": spec.dim, "tag": spec.tag}
    sel = {"name": str(spec["name"])}
    for key in ("kind", "dim", "tag"):
        sel[key] = spec.get(key)
    if sel["kind"] is not None and sel["kind"] not in KINDS:
        raise ValueError(f"{_PREFIX}unknown region kind '{sel['kind']}'")
    return sel


def _describe(r):
    return f"'{r.name}' ({r.kind}, dim {r.dim}, tag {r.tag})"


def _describe_sel(s):
    out = f"'{s['name']}'"
    for key in ("kind", "dim", "tag"):
        if s[key] is not None:
            out += f" {key} {s[key]}"
    return out


def _matches(r, s):
    return (
        r.name == s["name"]
        and (s["kind"] is None or r.kind == s["kind"])
        and (s["dim"] is None or r.dim == int(s["dim"]))
        and (s["tag"] is None or r.tag == int(s["tag"]))
    )


def _find(regions, spec):
    s = _selector(spec)
    hits = [i for i, r in enumerate(regions) if _matches(r, s)]
    if len(hits) == 1:
        return hits[0]
    if not hits:
        listed = ", ".join(_describe(r) for r in regions) or "none"
        raise ValueError(
            f"{_PREFIX}no region matches {_describe_sel(s)} (regions: {listed})"
        )
    listed = ", ".join(_describe(regions[i]) for i in hits)
    raise ValueError(
        f"{_PREFIX}{_describe_sel(s)} matches {len(hits)} regions ({listed}); "
        "pin the kind, dim or tag"
    )


def _rows(r):
    e = np.asarray(r.entries, dtype=np.int64)
    return e.reshape(-1, 2) if r.kind == "side" else e.reshape(-1, 1)


def _combine(op, first, rest):
    acc = {tuple(row) for row in _rows(first)}
    for r in rest:
        rows = {tuple(row) for row in _rows(r)}
        if op == "union":
            acc |= rows
        elif op == "intersection":
            acc &= rows
        else:
            acc -= rows
    width = 2 if first.kind == "side" else 1
    arr = np.array(sorted(acc), dtype=np.int64).reshape(-1, width)
    return arr if width == 2 else arr.reshape(-1)


def _check_free(regions, out, inputs, op):
    for i, r in enumerate(regions):
        if r.key == out.key and i not in inputs:
            raise ValueError(
                f"{_PREFIX}{op}: the result {_describe(out)} would replace an existing region"
            )


def _add(regions, region):
    for i, r in enumerate(regions):
        if r.key == region.key:
            regions[i] = region
            return
    regions.append(region)


def _apply(mesh, regions, edit):
    op = _op_name(edit.get("op"))
    inputs = [_find(regions, spec) for spec in edit.get("inputs", ())]
    nin = len(inputs)
    output = edit.get("output")
    out_dim = edit.get("dim")
    out_tag = edit.get("tag")

    if op == "delete":
        if nin == 0:
            raise ValueError(f"{_PREFIX}delete: no region given")
        return [r for i, r in enumerate(regions) if i not in set(inputs)]

    if op in ("rename", "retag"):
        if nin != 1:
            raise ValueError(f"{_PREFIX}{op}: takes exactly one region, got {nin}")
        if op == "rename" and not output:
            raise ValueError(f"{_PREFIX}rename: no new name given")
        if op == "retag" and out_tag is None and out_dim is None:
            raise ValueError(f"{_PREFIX}retag: give a new tag and/or dimension")
        src = regions[inputs[0]]
        out = Region(
            output if output else src.name,
            src.kind,
            src.entries.copy(),
            int(out_dim) if out_dim is not None else src.dim,
            int(out_tag) if out_tag is not None else src.tag,
        )
        _check_free(regions, out, inputs, op)
        if op == "retag" and out.kind == "cell" and "gmsh:physical" in mesh.cell_data:
            from ._common import warn

            warn(
                f"edit_regions: retagged cell region '{out.name}', but the mesh carries "
                "'gmsh:physical' cell data, which a gmsh write takes the physical tag "
                "from; that array still holds the old tag"
            )
        kept = [r for i, r in enumerate(regions) if i != inputs[0]]
        _add(kept, out)
        return kept

    if nin < 2:
        raise ValueError(f"{_PREFIX}{op}: takes two or more regions, got {nin}")
    if not output:
        raise ValueError(f"{_PREFIX}{op}: no name given for the result")
    first = regions[inputs[0]]
    rest = [regions[i] for i in inputs[1:]]
    dim = first.dim
    for r in rest:
        if r.kind != first.kind:
            raise ValueError(
                f"{_PREFIX}{op}: {_describe(first)} and {_describe(r)} are of different kinds"
            )
        if r.dim != dim:
            dim = -1
    out = Region(
        output,
        first.kind,
        _combine(op, first, rest),
        int(out_dim) if out_dim is not None else dim,
        int(out_tag) if out_tag is not None else -1,
    )
    _check_free(regions, out, inputs, op)
    kept = list(regions)
    if not edit.get("keep_inputs", True):
        kept = [r for i, r in enumerate(regions) if i not in set(inputs)]
    _add(kept, out)
    return kept


def edit_regions(mesh, edits):
    """Apply region edits, in order, to a copy of ``mesh``.

    :param mesh: the mesh; it is not modified.
    :param edits: a list of edit dicts -- ``op`` (``union``, ``intersection``,
        ``difference``, ``rename``, ``retag``, ``delete``), ``inputs`` (region
        names or ``{"name", "kind", "dim", "tag"}`` selectors), and as the
        operation needs ``output`` (the result's name), ``dim``, ``tag`` and
        ``keep_inputs`` (set operations; default ``True``). A single dict is
        accepted as a one-edit list.
    :returns: a copy of ``mesh`` with the edited regions, in the core's
        ``(kind, name, dim, tag)`` order.
    :raises ValueError: on a missing or ambiguous region, a set operation over
        different kinds, a missing name, or a result that would replace an
        unrelated region.
    """
    if isinstance(edits, dict):
        edits = [edits]
    out = mesh.copy()
    regions = [r.copy() for r in (getattr(mesh, "regions", None) or [])]
    for edit in edits:
        regions = _apply(mesh, regions, edit)
    regions.sort(key=lambda r: r.key)
    out.regions = regions
    return out


__all__ = ["edit_regions"]
