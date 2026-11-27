"""Helpers for complex-valued fields stored as ``<base>_real`` / ``<base>_imag`` pairs.

Public API:
    complex_pairs
    to_complex
    from_complex
"""

from __future__ import annotations

import numpy as np

REAL_SUFFIX = "_real"
IMAG_SUFFIX = "_imag"


def find_complex_pairs(names) -> list[tuple[str, str, str]]:
    """Return sorted ``(base, real_name, imag_name)`` tuples from ``names``."""
    present = set(names)
    pairs = []
    for name in names:
        if len(name) > len(REAL_SUFFIX) and name.endswith(REAL_SUFFIX):
            base = name[: -len(REAL_SUFFIX)]
            imag = base + IMAG_SUFFIX
            if imag in present:
                pairs.append((base, name, imag))
    pairs.sort(key=lambda item: item[0])
    return pairs


def complex_pairs(mesh, *, location: str = "both") -> dict[str, list[str]]:
    """List matched ``<base>_real`` / ``<base>_imag`` complex field base names.

    Args:
        mesh: A :class:`meshioplusplus.Mesh`.
        location: ``"point"``, ``"cell"``, or ``"both"`` (default).

    Returns:
        A dict with keys ``"point_data"`` and/or ``"cell_data"`` mapping to
        sorted lists of complex field base names.
    """
    if location not in ("point", "cell", "both"):
        raise ValueError(
            f"complex_pairs: location must be 'point', 'cell', or 'both', got {location!r}"
        )
    out: dict[str, list[str]] = {}
    if location in ("point", "both"):
        out["point_data"] = [
            base for base, _, _ in find_complex_pairs(mesh.point_data.keys())
        ]
    if location in ("cell", "both"):
        out["cell_data"] = [
            base for base, _, _ in find_complex_pairs(mesh.cell_data.keys())
        ]
    return out


def _combine_pair(re, im):
    re_arr = np.asarray(re)
    im_arr = np.asarray(im)
    if re_arr.shape != im_arr.shape:
        raise ValueError(
            f"Complex pair shape mismatch: real {re_arr.shape} vs imag {im_arr.shape}"
        )
    if re_arr.dtype == np.float32 and im_arr.dtype == np.float32:
        out = np.empty(re_arr.shape, dtype=np.complex64)
    else:
        out = np.empty(re_arr.shape, dtype=np.complex128)
    out.real = re_arr
    out.imag = im_arr
    return out


def to_complex(mesh, name: str | None = None, *, location: str = "both"):
    """Combine ``<name>_real`` and ``<name>_imag`` arrays into complex numpy arrays.

    Args:
        mesh: A :class:`meshioplusplus.Mesh`.
        name: Optional base name. When given, returns the complex ``np.ndarray``
            (for ``point_data``) or ``list[np.ndarray]`` (for ``cell_data``).
            When ``None``, returns a dict of all complex pairs in ``location``.
        location: ``"point"``, ``"cell"``, or ``"both"`` (default).

    Returns:
        Complex numpy array(s) with dtype ``complex64`` when both inputs are
        ``float32``, and ``complex128`` otherwise.
    """
    if location not in ("point", "cell", "both"):
        raise ValueError(
            f"to_complex: location must be 'point', 'cell', or 'both', got {location!r}"
        )
    if name is not None:
        re_name = name + REAL_SUFFIX
        im_name = name + IMAG_SUFFIX
        if location in ("point", "both"):
            if re_name in mesh.point_data and im_name in mesh.point_data:
                return _combine_pair(mesh.point_data[re_name], mesh.point_data[im_name])
        if location in ("cell", "both"):
            if re_name in mesh.cell_data and im_name in mesh.cell_data:
                re_blocks = mesh.cell_data[re_name]
                im_blocks = mesh.cell_data[im_name]
                if len(re_blocks) != len(im_blocks):
                    raise ValueError(
                        f"Complex cell_data pair '{name}' block count mismatch: "
                        f"{len(re_blocks)} vs {len(im_blocks)}"
                    )
                return [_combine_pair(r, i) for r, i in zip(re_blocks, im_blocks)]
        raise KeyError(f"No complex pair '{re_name}' / '{im_name}' found in {location}")

    out: dict[str, dict] = {}
    if location in ("point", "both"):
        out["point_data"] = {
            base: _combine_pair(mesh.point_data[re_k], mesh.point_data[im_k])
            for base, re_k, im_k in find_complex_pairs(mesh.point_data.keys())
        }
    if location in ("cell", "both"):
        out["cell_data"] = {
            base: [
                _combine_pair(r, i)
                for r, i in zip(mesh.cell_data[re_k], mesh.cell_data[im_k])
            ]
            for base, re_k, im_k in find_complex_pairs(mesh.cell_data.keys())
        }
    return out


def from_complex(mesh, name: str, array, *, location: str = "point"):
    """Store a complex array on ``mesh`` as ``<name>_real`` and ``<name>_imag``.

    Args:
        mesh: A :class:`meshioplusplus.Mesh` (modified in place and returned).
        name: Base field name (without ``_real`` / ``_imag``).
        array: A complex ``np.ndarray`` (for ``location="point"``) or a list of
            complex ``np.ndarray``s (for ``location="cell"``).
        location: ``"point"`` (default) or ``"cell"``.

    Returns:
        ``mesh``, for chaining.
    """
    if not name:
        raise ValueError("from_complex: name must be non-empty")
    re_name = name + REAL_SUFFIX
    im_name = name + IMAG_SUFFIX
    if location == "point":
        arr = np.asarray(array)
        mesh.point_data[re_name] = np.ascontiguousarray(np.real(arr))
        mesh.point_data[im_name] = np.ascontiguousarray(np.imag(arr))
    elif location == "cell":
        blocks = [np.asarray(b) for b in array]
        mesh.cell_data[re_name] = [np.ascontiguousarray(np.real(b)) for b in blocks]
        mesh.cell_data[im_name] = [np.ascontiguousarray(np.imag(b)) for b in blocks]
    else:
        raise ValueError(
            f"from_complex: location must be 'point' or 'cell', got {location!r}"
        )
    return mesh
