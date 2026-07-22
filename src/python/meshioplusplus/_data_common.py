"""Shared vocabulary for the *data* operations.

These operations act on a mesh's ``point_data`` / ``cell_data`` / ``field_data``
arrays rather than on its geometry, which they never modify. A dependency-free
mesh *operation* (not a file format), so nothing here is registered in the
format registry.

Non-finite values (NaN, +/-inf) follow one policy across every data operation:
they are **always excluded from reductions** (min/max/mean/stddev and the
averaging accumulators), and ``nan_policy`` only decides what reaches the
output -- ``"ignore"`` passes them through, ``"replace"`` substitutes
``nan_replacement``, and ``"fail"`` raises. ``data_info`` never raises; it
counts them instead.

Public API:
    normalize_location -- canonicalize a location name
    location_map       -- the mesh attribute holding a location's arrays
    available_keys     -- sorted names at a location
    require_key        -- raise a listing error if a key is absent
"""

from __future__ import annotations

_LOCATION_ALIASES = {
    "point": "point_data",
    "point_data": "point_data",
    "node": "point_data",
    "cell": "cell_data",
    "cell_data": "cell_data",
    "element": "cell_data",
    "field": "field_data",
    "field_data": "field_data",
}

LOCATIONS = ("point_data", "cell_data", "field_data")


def normalize_location(location: str) -> str:
    """Return the canonical ``*_data`` name for ``location``."""
    try:
        return _LOCATION_ALIASES[location]
    except KeyError:
        raise ValueError(
            f"meshio++: unknown data location '{location}' "
            "(expected 'point', 'cell' or 'field')"
        ) from None


def location_map(mesh, location: str) -> dict:
    """Return the mesh dict holding ``location``'s arrays."""
    return getattr(mesh, normalize_location(location))


def available_keys(mesh, location: str) -> list:
    """Return the sorted array names present at ``location``."""
    return sorted(location_map(mesh, location).keys())


def require_key(mesh, location: str, name: str) -> None:
    """Raise ``ValueError`` naming the available keys if ``name`` is absent.

    ``ValueError`` rather than the more idiomatic ``KeyError`` so that both code
    paths agree: the C++ core raises ``std::invalid_argument``, which pybind11
    surfaces as ``ValueError``. Windows CI runs entirely on the numpy fallbacks,
    so a mismatch here would be a platform-dependent test failure.
    """
    loc = normalize_location(location)
    if name in location_map(mesh, loc):
        return
    keys = available_keys(mesh, loc)
    detail = f"available: {', '.join(keys)}" if keys else f"the mesh has no {loc}"
    raise ValueError(f"meshio++: no {loc} array named '{name}' ({detail})")


def num_components(array) -> int:
    """Return the product of ``array``'s trailing dimensions (1 for a scalar)."""
    import numpy as np

    a = np.asarray(array)
    if a.ndim < 2:
        return 1
    n = 1
    for d in a.shape[1:]:
        n *= int(d)
    return max(n, 1)
