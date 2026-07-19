"""Mesh comparison ("diff").

A dependency-free mesh *operation* (not a file format): compare two meshes and
report whether they are equivalent within a tolerance, plus a structured
description of any differences. The C++ core (``_core.diff``) does the work;
this module is the thin shim (try C++, fall back to a pure-numpy reference)
plus the ``point_sets``/``cell_sets`` comparison the core cannot see across the
Python boundary, and a compact human-readable formatter for the CLI.

The elementwise tolerance rule is ``abs_err <= atol + rtol*|expected|`` with
``expected`` taken from the first mesh. The overall ``verdict`` is one of
``"identical"`` (structurally equal and every value bitwise-equal),
``"equal within tolerance"`` (values drifted but within tolerance), or
``"different"``.

Public API:

* :func:`diff` — full structured report (a nested dict).
* :func:`meshes_equal` — convenience boolean, handy in test suites.
"""

from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------- #
# elementwise array comparison                                                #
# --------------------------------------------------------------------------- #
def _empty_array_diff(name=""):
    return {
        "name": name,
        "shape_mismatch": False,
        "size_a": 0,
        "size_b": 0,
        "max_abs_error": 0.0,
        "max_rel_error": 0.0,
        "worst_index": -1,
        "num_exceeding": 0,
        "exact": True,
    }


def _compare_array(a, b, row_map, atol, rtol, name=""):
    """Row-wise compare two arrays; ``row_map`` maps B rows -> A rows (or None)."""
    a = np.asarray(a)
    b = np.asarray(b)
    out = _empty_array_diff(name)
    out["size_a"] = int(a.size)
    out["size_b"] = int(b.size)

    a2 = a.reshape(a.shape[0], -1) if a.ndim >= 1 and a.shape[0] else a.reshape(1, -1)
    b2 = b.reshape(b.shape[0], -1) if b.ndim >= 1 and b.shape[0] else b.reshape(1, -1)
    rows_a, cols_a = a2.shape if a.size else (0, a2.shape[1] if a2.ndim == 2 else 1)
    rows_b, cols_b = b2.shape if b.size else (0, b2.shape[1] if b2.ndim == 2 else 1)

    if cols_a != cols_b or (row_map is None and rows_a != rows_b):
        out["shape_mismatch"] = True
        out["exact"] = False
        return out
    nrows = rows_b if row_map is not None else rows_a
    ncols = cols_a
    if nrows == 0 or ncols == 0:
        return out

    exp = a2[np.asarray(row_map, dtype=np.int64)] if row_map is not None else a2[:nrows]
    obs = b2[:nrows]
    exp = exp.astype(np.float64, copy=False)
    obs = obs.astype(np.float64, copy=False)
    absd = np.abs(exp - obs)
    denom = np.abs(exp)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(denom > 0.0, absd / denom, 0.0)
    exceed = absd > (atol + rtol * denom)

    flat = absd.reshape(-1)
    worst = int(np.argmax(flat)) if flat.size else -1
    out["max_abs_error"] = float(flat[worst]) if worst >= 0 else 0.0
    out["max_rel_error"] = float(rel.max()) if rel.size else 0.0
    out["worst_index"] = worst if flat.size else -1
    out["num_exceeding"] = int(exceed.sum())
    out["exact"] = bool(np.array_equal(exp, obs))
    return out


# --------------------------------------------------------------------------- #
# cell block comparison                                                       #
# --------------------------------------------------------------------------- #
def _cell_node_lists(blk, node_map):
    """Yield each cell's node-id list, optionally remapped through ``node_map``."""
    data = blk.data

    def remap(ids):
        if node_map is None:
            return [int(x) for x in ids]
        n = len(node_map)
        return [int(node_map[int(x)]) if 0 <= int(x) < n else int(x) for x in ids]

    if isinstance(data, np.ndarray) and data.ndim == 2:
        for row in data:
            yield remap(row)
    else:
        for cell in data:
            if (
                isinstance(cell, list)
                and cell
                and isinstance(cell[0], (list, tuple, np.ndarray))
            ):
                flat = []
                for face in cell:
                    flat.extend(int(x) for x in np.asarray(face).reshape(-1))
                yield remap(flat)
            else:
                yield remap(np.asarray(cell).reshape(-1))


def _compare_block(blk_a, blk_b, block, node_map, max_reported):
    out = {
        "block": int(block),
        "type_a": blk_a.type,
        "type_b": blk_b.type,
        "count_a": int(len(blk_a.data)),
        "count_b": int(len(blk_b.data)),
        "type_mismatch": blk_a.type != blk_b.type,
        "count_mismatch": len(blk_a.data) != len(blk_b.data),
        "conn_mismatch_count": 0,
        "first_mismatches": [],
    }
    if out["type_mismatch"] or out["count_mismatch"]:
        return out
    la = _cell_node_lists(blk_a, None)
    lb = _cell_node_lists(blk_b, node_map)
    for c, (na, nb) in enumerate(zip(la, lb)):
        if list(na) != list(nb):
            out["conn_mismatch_count"] += 1
            if len(out["first_mismatches"]) < max_reported:
                out["first_mismatches"].append(int(c))
    return out


# --------------------------------------------------------------------------- #
# data section comparison                                                     #
# --------------------------------------------------------------------------- #
def _compare_point_data(a, b, row_map, atol, rtol):
    names_a = sorted(a.point_data.keys())
    names_b = sorted(b.point_data.keys())
    sa, sb = set(names_a), set(names_b)
    out = {
        "only_in_a": sorted(sa - sb),
        "only_in_b": sorted(sb - sa),
        "shared": [],
    }
    for name in sorted(sa & sb):
        out["shared"].append(
            _compare_array(
                a.point_data[name], b.point_data[name], row_map, atol, rtol, name
            )
        )
    return out


def _fold(into, add, base):
    into["size_a"] += add["size_a"]
    into["size_b"] += add["size_b"]
    into["shape_mismatch"] = into["shape_mismatch"] or add["shape_mismatch"]
    into["num_exceeding"] += add["num_exceeding"]
    into["exact"] = into["exact"] and add["exact"]
    into["max_rel_error"] = max(into["max_rel_error"], add["max_rel_error"])
    if add["worst_index"] != -1 and add["max_abs_error"] > into["max_abs_error"]:
        into["max_abs_error"] = add["max_abs_error"]
        into["worst_index"] = base + add["worst_index"]


def _compare_cell_data(a, b, atol, rtol):
    names_a = sorted(a.cell_data.keys())
    names_b = sorted(b.cell_data.keys())
    sa, sb = set(names_a), set(names_b)
    out = {
        "only_in_a": sorted(sa - sb),
        "only_in_b": sorted(sb - sa),
        "shared": [],
    }
    for name in sorted(sa & sb):
        blocks_a = a.cell_data[name]
        blocks_b = b.cell_data[name]
        summary = _empty_array_diff(name)
        if len(blocks_a) != len(blocks_b):
            summary["shape_mismatch"] = True
            summary["exact"] = False
        else:
            base = 0
            for ba, bb in zip(blocks_a, blocks_b):
                if ba is None or bb is None:
                    continue
                part = _compare_array(ba, bb, None, atol, rtol, name)
                _fold(summary, part, base)
                base += int(np.asarray(ba).size)
        out["shared"].append(summary)
    return out


def _compare_field_data(a, b, atol, rtol):
    names_a = sorted(a.field_data.keys())
    names_b = sorted(b.field_data.keys())
    sa, sb = set(names_a), set(names_b)
    out = {
        "only_in_a": sorted(sa - sb),
        "only_in_b": sorted(sb - sa),
        "shared": [],
    }
    for name in sorted(sa & sb):
        out["shared"].append(
            _compare_array(
                a.field_data[name], b.field_data[name], None, atol, rtol, name
            )
        )
    return out


# --------------------------------------------------------------------------- #
# unordered node matching (bucket grid, no O(N^2) scan)                       #
# --------------------------------------------------------------------------- #
def _match_points(a, b, atol, rtol):
    """Match each B point to a unique A point by proximity. Returns (a_of_b,
    max_residual, reason) — a_of_b is None on failure."""
    pa = np.asarray(a.points, dtype=np.float64)
    pb = np.asarray(b.points, dtype=np.float64)
    if pa.ndim == 1:
        pa = pa.reshape(len(pa), -1)
    if pb.ndim == 1:
        pb = pb.reshape(len(pb), -1)
    na, nb = pa.shape[0], pb.shape[0]
    if na != nb:
        return None, 0.0, f"point counts differ ({na} vs {nb})"
    if pa.shape[1] != pb.shape[1]:
        return None, 0.0, "point dimensions differ"
    pdim = pa.shape[1]
    gdim = min(pdim, 3)

    lo = pa[:, :gdim].min(axis=0) if na else np.zeros(gdim)
    hi = pa[:, :gdim].max(axis=0) if na else np.zeros(gdim)
    h = np.empty(gdim)
    for d in range(gdim):
        cell = atol + rtol * max(abs(lo[d]), abs(hi[d]))
        if not cell > 0.0:
            cell = (hi[d] - lo[d]) if hi[d] > lo[d] else 1.0
        h[d] = cell

    def bucket(p):
        return tuple(
            int(np.floor((p[d] - lo[d]) / h[d])) if d < gdim else 0 for d in range(3)
        )

    grid = {}
    for i in range(na):
        grid.setdefault(bucket(pa[i]), []).append(i)

    a_of_b = np.full(nb, -1, dtype=np.int64)
    used = bytearray(na)
    max_res = 0.0
    for j in range(nb):
        bj = bucket(pb[j])
        best, best_err = -1, np.inf
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for ai in grid.get((bj[0] + dx, bj[1] + dy, bj[2] + dz), ()):
                        if used[ai]:
                            continue
                        e = np.abs(pa[ai] - pb[j])
                        if np.all(e <= atol + rtol * np.abs(pa[ai])):
                            werr = float(e.max())
                            if werr < best_err:
                                best_err, best = werr, ai
        if best < 0:
            return (
                None,
                0.0,
                (
                    f"point {j} of the second mesh has no unmatched neighbour within tolerance"
                ),
            )
        used[best] = 1
        a_of_b[j] = best
        max_res = max(max_res, best_err)
    return a_of_b, max_res, ""


# --------------------------------------------------------------------------- #
# verdict                                                                     #
# --------------------------------------------------------------------------- #
_VERDICTS = ("identical", "equal within tolerance", "different")


def _array_verdict(ad):
    if ad["shape_mismatch"] or ad["num_exceeding"] > 0:
        return 2
    if not ad["exact"]:
        return 1
    return 0


def _recompute_verdict(report):
    v = 0
    if report["point_count_mismatch"] or report["block_count_mismatch"]:
        v = 2
    v = max(v, _array_verdict(report["points"]))
    for bd in report["blocks"]:
        if bd["type_mismatch"] or bd["count_mismatch"] or bd["conn_mismatch_count"] > 0:
            v = 2
    for section in ("point_data", "cell_data", "field_data"):
        dd = report[section]
        if dd["only_in_a"] or dd["only_in_b"]:
            v = 2
        for ad in dd["shared"]:
            v = max(v, _array_verdict(ad))
    return _VERDICTS[v]


# --------------------------------------------------------------------------- #
# pure-numpy reference (fallback when the C++ core is unavailable)            #
# --------------------------------------------------------------------------- #
def _diff_py(a, b, atol, rtol, unordered, max_reported):
    report = {
        "verdict": "identical",
        "unordered": unordered,
        "correspondence_failed": False,
        "point_count_mismatch": False,
        "num_points_a": int(len(a.points)),
        "num_points_b": int(len(b.points)),
        "points": _empty_array_diff(),
        "block_count_mismatch": len(a.cells) != len(b.cells),
        "num_blocks_a": int(len(a.cells)),
        "num_blocks_b": int(len(b.cells)),
        "blocks": [],
        "point_data": {"only_in_a": [], "only_in_b": [], "shared": []},
        "cell_data": {"only_in_a": [], "only_in_b": [], "shared": []},
        "field_data": {"only_in_a": [], "only_in_b": [], "shared": []},
        "messages": [],
    }

    row_map = None
    node_map = None
    if unordered:
        a_of_b, residual, reason = _match_points(a, b, atol, rtol)
        if a_of_b is None:
            report["correspondence_failed"] = True
            report["verdict"] = "different"
            report["messages"].append(f"unordered matching failed: {reason}")
            report["point_count_mismatch"] = len(a.points) != len(b.points)
            return report
        row_map = a_of_b
        node_map = a_of_b
        report["points"]["max_abs_error"] = residual
        report["points"]["exact"] = residual == 0.0
        report["messages"].append(
            f"unordered: matched {report['num_points_b']} points by proximity"
        )

    pa = np.asarray(a.points)
    pb = np.asarray(b.points)
    dim_a = pa.shape[1] if pa.ndim == 2 else 1
    dim_b = pb.shape[1] if pb.ndim == 2 else 1
    if len(a.points) != len(b.points) or dim_a != dim_b:
        report["point_count_mismatch"] = True
        report["points"]["shape_mismatch"] = True
        report["points"]["exact"] = False
        report["points"]["size_a"] = int(pa.size)
        report["points"]["size_b"] = int(pb.size)
    else:
        report["points"] = _compare_array(pa, pb, row_map, atol, rtol)

    nblocks = min(len(a.cells), len(b.cells))
    for i in range(nblocks):
        report["blocks"].append(
            _compare_block(a.cells[i], b.cells[i], i, node_map, max_reported)
        )

    report["point_data"] = _compare_point_data(a, b, row_map, atol, rtol)
    report["cell_data"] = _compare_cell_data(a, b, atol, rtol)
    report["field_data"] = _compare_field_data(a, b, atol, rtol)
    report["verdict"] = _recompute_verdict(report)
    return report


# --------------------------------------------------------------------------- #
# sets comparison (Python-side; the C++ core cannot see point/cell sets)      #
# --------------------------------------------------------------------------- #
def _compare_sets(a, b):
    ps_a = getattr(a, "point_sets", None) or {}
    ps_b = getattr(b, "point_sets", None) or {}
    cs_a = getattr(a, "cell_sets", None) or {}
    cs_b = getattr(b, "cell_sets", None) or {}

    def point_changed(name):
        return set(np.asarray(ps_a[name]).tolist()) != set(
            np.asarray(ps_b[name]).tolist()
        )

    def cell_changed(name):
        ba, bb = cs_a[name], cs_b[name]
        if len(ba) != len(bb):
            return True
        for xa, xb in zip(ba, bb):
            xa = np.asarray([]) if xa is None else np.asarray(xa)
            xb = np.asarray([]) if xb is None else np.asarray(xb)
            if set(xa.tolist()) != set(xb.tolist()):
                return True
        return False

    return {
        "point_sets": {
            "only_in_a": sorted(set(ps_a) - set(ps_b)),
            "only_in_b": sorted(set(ps_b) - set(ps_a)),
            "changed": sorted(n for n in (set(ps_a) & set(ps_b)) if point_changed(n)),
        },
        "cell_sets": {
            "only_in_a": sorted(set(cs_a) - set(cs_b)),
            "only_in_b": sorted(set(cs_b) - set(cs_a)),
            "changed": sorted(n for n in (set(cs_a) & set(cs_b)) if cell_changed(n)),
        },
    }


def _sets_differ(sets):
    for section in ("point_sets", "cell_sets"):
        s = sets[section]
        if s["only_in_a"] or s["only_in_b"] or s["changed"]:
            return True
    return False


# --------------------------------------------------------------------------- #
# public API                                                                  #
# --------------------------------------------------------------------------- #
def diff(
    mesh_a, mesh_b, *, atol=1e-12, rtol=1e-9, unordered=False, max_reported=10
) -> dict:
    """Compare two meshes and return a structured report.

    Parameters
    ----------
    mesh_a, mesh_b :
        the meshes to compare (``mesh_a`` provides the "expected" values).
    atol, rtol :
        tolerances in ``abs_err <= atol + rtol*|expected|`` (defaults
        ``1e-12`` / ``1e-9``).
    unordered :
        match points by spatial proximity instead of by index (tolerant to a
        shuffled node order). Conservative: reports ``correspondence_failed``
        rather than guessing when a one-to-one match can't be built.
    max_reported :
        cap on first-mismatching cell indices recorded per block.

    Returns
    -------
    dict
        A nested report with an overall ``"verdict"`` (``"identical"``,
        ``"equal within tolerance"``, or ``"different"``), a ``"points"``
        summary, a ``"blocks"`` list, ``"point_data"``/``"cell_data"``/
        ``"field_data"`` sections, a ``"sets"`` section (Python-only), and
        ``"messages"``.
    """
    report = None
    try:
        from . import _core

        report = _core.diff(mesh_a, mesh_b, atol, rtol, unordered, int(max_reported))
    except Exception:
        report = None
    if report is None:
        report = _diff_py(mesh_a, mesh_b, atol, rtol, unordered, int(max_reported))

    # Named point_sets/cell_sets live only on the Python Mesh; compare them here
    # and fold any difference into the verdict.
    sets = _compare_sets(mesh_a, mesh_b)
    report["sets"] = sets
    if _sets_differ(sets):
        if report["verdict"] != "different":
            report.setdefault("messages", []).append(
                "named point_sets/cell_sets differ"
            )
        report["verdict"] = "different"
    return report


def meshes_equal(mesh_a, mesh_b, *, atol=1e-12, rtol=1e-9, unordered=False) -> bool:
    """Are two meshes equal within tolerance? (Geometry + data; not sets.)

    A convenient boolean for test suites — ``True`` when the verdict would be
    ``"identical"`` or ``"equal within tolerance"``. Named ``point_sets`` /
    ``cell_sets`` are ignored (use :func:`diff` to compare those).
    """
    if not unordered:
        try:
            from . import _core

            return bool(_core.meshes_equal(mesh_a, mesh_b, atol, rtol))
        except Exception:
            pass
    report = diff(mesh_a, mesh_b, atol=atol, rtol=rtol, unordered=unordered)
    # diff() may flip the verdict to "different" purely on sets; recompute a
    # sets-agnostic verdict for the boolean helper.
    verdict = (
        _recompute_verdict(report)
        if not report["correspondence_failed"]
        else "different"
    )
    return verdict != "different"


# --------------------------------------------------------------------------- #
# human-readable formatting (for the CLI)                                     #
# --------------------------------------------------------------------------- #
def _format_report(report) -> str:
    """Compact, differences-first rendering of a diff report."""
    lines = []
    lines.append(f"verdict: {report['verdict']}")
    for msg in report.get("messages", []):
        lines.append(f"  note: {msg}")

    if report["point_count_mismatch"]:
        lines.append(
            f"points: count/dim differ ({report['num_points_a']} vs {report['num_points_b']})"
        )
    else:
        p = report["points"]
        if p["shape_mismatch"]:
            lines.append("points: shape mismatch")
        elif not p["exact"]:
            lines.append(
                f"points: max_abs={p['max_abs_error']:.3e} max_rel={p['max_rel_error']:.3e} "
                f"exceeding={p['num_exceeding']} worst_index={p['worst_index']}"
            )

    if report["block_count_mismatch"]:
        lines.append(
            f"cells: block count differs ({report['num_blocks_a']} vs {report['num_blocks_b']})"
        )
    for bd in report["blocks"]:
        if bd["type_mismatch"]:
            lines.append(f"cells[{bd['block']}]: type {bd['type_a']} vs {bd['type_b']}")
        elif bd["count_mismatch"]:
            lines.append(
                f"cells[{bd['block']}] ({bd['type_a']}): count {bd['count_a']} vs {bd['count_b']}"
            )
        elif bd["conn_mismatch_count"] > 0:
            first = ", ".join(str(x) for x in bd["first_mismatches"])
            lines.append(
                f"cells[{bd['block']}] ({bd['type_a']}): "
                f"{bd['conn_mismatch_count']} connectivity mismatch(es) [first: {first}]"
            )

    for section in ("point_data", "cell_data", "field_data"):
        dd = report[section]
        if dd["only_in_a"]:
            lines.append(f"{section}: only in A: {', '.join(dd['only_in_a'])}")
        if dd["only_in_b"]:
            lines.append(f"{section}: only in B: {', '.join(dd['only_in_b'])}")
        for ad in dd["shared"]:
            if ad["shape_mismatch"]:
                lines.append(f"{section}['{ad['name']}']: shape mismatch")
            elif not ad["exact"]:
                lines.append(
                    f"{section}['{ad['name']}']: max_abs={ad['max_abs_error']:.3e} "
                    f"max_rel={ad['max_rel_error']:.3e} exceeding={ad['num_exceeding']}"
                )

    sets = report.get("sets")
    if sets:
        for section in ("point_sets", "cell_sets"):
            s = sets[section]
            if s["only_in_a"]:
                lines.append(f"{section}: only in A: {', '.join(s['only_in_a'])}")
            if s["only_in_b"]:
                lines.append(f"{section}: only in B: {', '.join(s['only_in_b'])}")
            if s["changed"]:
                lines.append(f"{section}: changed: {', '.join(s['changed'])}")

    return "\n".join(lines)
