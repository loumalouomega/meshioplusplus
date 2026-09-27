"""Pass/fail thresholds over ``compute_quality``'s metrics -- a quality gate.

``check_quality(mesh, "scaled_jacobian >= 0.2; aspect_ratio <= 5 @ 1%")``
scores every cell and tests each threshold against the per-cell values (the
report's histograms span each metric's own range, so they cannot say how many
cells fall below a bound). Cells where a metric does not apply (NaN) are not
evaluated; inverted and degenerate cells are gated by count. The semantics are
``operations/quality_gate.hpp``'s; the numpy twin below runs on the Python
``compute_quality`` when the core is unavailable.

The specification text: clauses separated by ``;`` or ``,``, each
``METRIC >= VALUE`` or ``METRIC <= VALUE`` (a newline separates clauses too,
and ``#`` starts a comment, so a gate file is one clause per line), optionally
followed by
``@ FRACTION`` (a number in ``[0, 1]`` or a percentage such as ``@1%``) --
the fraction of evaluated cells allowed to violate the bound.
"""

from __future__ import annotations

import math
import re

import numpy as np

from ._fallback import core_op_declined

_PREFIX = "meshio++: quality gate: "

#: ``compute_quality``'s metrics, in its fixed order.
METRICS = (
    "quality:volume",
    "quality:scaled_jacobian",
    "quality:aspect_ratio",
    "quality:skewness",
    "quality:min_angle",
    "quality:max_angle",
    "quality:warpage",
    "quality:min_dihedral",
    "quality:max_dihedral",
    "quality:inverted",
    "quality:degenerate",
)


def _full_name(metric):
    full = metric if metric.startswith("quality:") else "quality:" + metric
    if full not in METRICS:
        names = ", ".join(m[8:] for m in METRICS)
        raise ValueError(
            f"{_PREFIX}unknown metric '{metric}' (expected one of {names})"
        )
    return full


def _number(text, clause):
    text = text.strip()
    try:
        value = float(text)
    except ValueError:
        value = math.nan
    if not text or not math.isfinite(value):
        raise ValueError(f"{_PREFIX}'{clause}': expected a number, got '{text}'")
    return value


def parse_quality_thresholds(spec: str) -> list:
    """The thresholds a specification text names, as dicts
    ``{"metric", "min", "max", "max_fraction"}``."""
    out = []
    spec = "\n".join(line.split("#", 1)[0] for line in spec.splitlines())
    for clause in re.split(r"[;,\n]", spec):
        clause = clause.strip()
        if not clause:
            continue
        ge, le = clause.find(">="), clause.find("<=")
        if (ge < 0) == (le < 0):
            raise ValueError(
                f"{_PREFIX}'{clause}': expected METRIC >= VALUE or METRIC <= VALUE"
            )
        op = ge if ge >= 0 else le
        t = {
            "metric": _full_name(clause[:op].strip()),
            "min": None,
            "max": None,
            "max_fraction": 0.0,
        }
        rest = clause[op + 2 :]
        if "@" in rest:
            rest, frac = rest.split("@", 1)
            frac = frac.strip()
            percent = frac.endswith("%")
            t["max_fraction"] = _number(frac.rstrip("%"), clause) / (
                100.0 if percent else 1.0
            )
        t["min" if ge >= 0 else "max"] = _number(rest, clause)
        out.append(t)
    return out


def _g(v):
    return f"{v:g}"


def _name(metric, lo, hi, frac):
    m = metric[8:]
    if lo is not None and hi is not None:
        s = f"{_g(lo)} <= {m} <= {_g(hi)}"
    elif lo is not None:
        s = f"{m} >= {_g(lo)}"
    else:
        s = f"{m} <= {_g(hi)}"
    if frac > 0.0:
        s += f" @ {_g(frac * 100.0)}%"
    return s


def _normalise(require):
    if require is None:
        return []
    if isinstance(require, str):
        return parse_quality_thresholds(require)
    out = []
    for item in require:
        if isinstance(item, str):
            out.extend(parse_quality_thresholds(item))
        else:
            out.append(
                {
                    "metric": item["metric"],
                    "min": item.get("min"),
                    "max": item.get("max"),
                    "max_fraction": float(item.get("max_fraction", 0.0) or 0.0),
                }
            )
    return out


def _check_quality_py(mesh, thresholds, max_inverted, max_degenerate):
    from ._quality import compute_quality

    rep = compute_quality(mesh)
    checks = []
    for t in thresholds:
        metric = _full_name(t["metric"])
        lo = t["min"] if t["min"] is not None and not math.isnan(t["min"]) else None
        hi = t["max"] if t["max"] is not None and not math.isnan(t["max"]) else None
        frac = float(t["max_fraction"])
        if lo is None and hi is None:
            raise ValueError(
                f"{_PREFIX}threshold on '{t['metric']}' has neither a minimum nor a maximum"
            )
        if not 0.0 <= frac <= 1.0:
            raise ValueError(
                f"{_PREFIX}threshold on '{t['metric']}': the allowed fraction must lie in [0, 1]"
            )
        arrays = rep["cell_arrays"].get(metric, [])
        values = (
            np.concatenate([np.asarray(a, dtype=float).reshape(-1) for a in arrays])
            if arrays
            else np.empty(0)
        )
        finite = np.isfinite(values)
        margin = np.full(len(values), np.inf)
        if lo is not None:
            margin = np.minimum(margin, values - lo)
        if hi is not None:
            margin = np.minimum(margin, hi - values)
        margin[~finite] = np.inf
        evaluated = int(finite.sum())
        violations = int((margin < 0.0).sum())
        worst, worst_cell = math.nan, -1
        if evaluated:
            k = int(np.argmin(np.where(finite, margin, np.inf)))  # first minimum
            worst, worst_cell = float(values[k]), k
        fraction = violations / evaluated if evaluated else 0.0
        if not evaluated:
            from ._common import warn

            warn(
                f"quality gate: '{_name(metric, lo, hi, frac)}' applies to no cell of this "
                "mesh (the metric is not defined for its cell types), so it passes vacuously"
            )
        checks.append(
            {
                "name": _name(metric, lo, hi, frac),
                "metric": metric,
                "min": math.nan if lo is None else float(lo),
                "max": math.nan if hi is None else float(hi),
                "max_fraction": frac,
                "evaluated": evaluated,
                "violations": violations,
                "fraction": fraction,
                "worst": worst,
                "worst_cell": worst_cell,
                "passed": violations == 0 or fraction <= frac,
            }
        )
    for what, count, limit in (
        ("inverted", rep["num_inverted"], max_inverted),
        ("degenerate", rep["num_degenerate"], max_degenerate),
    ):
        if limit < 0:
            continue
        n = rep["num_cells"]
        checks.append(
            {
                "name": f"{what} <= {int(limit)}",
                "metric": what,
                "min": math.nan,
                "max": float(limit),
                "max_fraction": 0.0,
                "evaluated": int(n),
                "violations": int(count),
                "fraction": count / n if n else 0.0,
                "worst": float(count),
                "worst_cell": -1,
                "passed": count <= limit,
            }
        )
    return {
        "passed": all(c["passed"] for c in checks),
        "num_cells": int(rep["num_cells"]),
        "num_inverted": int(rep["num_inverted"]),
        "num_degenerate": int(rep["num_degenerate"]),
        "checks": checks,
    }


def check_quality(
    mesh, require=None, max_inverted: int = 0, max_degenerate: int = 0
) -> dict:
    """Test a mesh's cell quality against thresholds.

    :param mesh: the mesh to score.
    :param require: the thresholds -- a specification text
        (``"scaled_jacobian >= 0.2; aspect_ratio <= 5 @ 1%"``), a list of such
        texts, or a list of ``{"metric", "min", "max", "max_fraction"}`` dicts.
        A metric is a ``compute_quality`` name, with or without ``quality:``.
    :param max_inverted: the most inverted cells allowed; negative disables it.
    :param max_degenerate: the most degenerate cells allowed; negative
        disables it.
    :returns: ``{"passed", "num_cells", "num_inverted", "num_degenerate",
        "checks"}``, each check a dict with ``name``, ``metric``, ``min``,
        ``max``, ``max_fraction``, ``evaluated``, ``violations``, ``fraction``,
        ``worst`` (the value closest to or furthest beyond the bound),
        ``worst_cell`` (global block-major index) and ``passed``.
    """
    thresholds = _normalise(require)
    try:
        from . import _core

        return dict(
            _core.check_quality(
                mesh, thresholds, int(max_inverted), int(max_degenerate)
            )
        )
    except Exception as exc:
        if not core_op_declined(exc, "check_quality"):
            raise
    return _check_quality_py(mesh, thresholds, int(max_inverted), int(max_degenerate))


def format_quality_gate(result: dict) -> str:
    """The human-readable summary both CLIs print (``quality_gate_summary``)."""
    lines = [
        f"quality gate: {'PASS' if result['passed'] else 'FAIL'} ({result['num_cells']} cells)"
    ]
    for c in result["checks"]:
        line = (
            f"  {'pass' if c['passed'] else 'FAIL'}  {c['name']}: {c['violations']} of "
            f"{c['evaluated']} cells violate"
        )
        if c["metric"].startswith("quality:") and c["worst_cell"] >= 0:
            line += f", worst {_g(c['worst'])} at cell {c['worst_cell']}"
        lines.append(line)
    return "\n".join(lines) + "\n"


__all__ = ["check_quality"]
