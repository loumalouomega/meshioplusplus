"""Array management for a mesh's data: rename, drop, and keep-only.

A dependency-free mesh *operation* (not a file format): it rewrites *which*
data arrays a mesh carries and under what names, copying values verbatim and
leaving the geometry (points and every cell block) bit-identical.

The three phases run in a fixed order: **keep**, then **drop**, then
**rename**. ``keep`` is a whitelist that only applies to the locations it
actually mentions, so keeping ``point_data`` arrays leaves ``cell_data`` and
``field_data`` untouched.

``point_sets`` / ``cell_sets`` pass through unchanged -- they index geometry,
which these operations never modify.

Public API:
    data_manage -- combined keep/drop/rename
    data_drop   -- drop the named arrays at one location
    data_keep   -- keep only the named arrays at one location
    data_rename -- rename one array
"""

from __future__ import annotations

import copy

from ._data_common import LOCATIONS, location_map, normalize_location, require_key


def _clone(mesh):
    """A deep copy of ``mesh`` that shares nothing with the original."""
    return copy.deepcopy(mesh)


def _manage_py(mesh, keep, drop, rename, ignore_missing):
    """Pure-Python reference implementation of :func:`data_manage`."""
    # Validate everything against the input before touching anything, so a bad
    # key costs no work -- the same contract the C++ core offers.
    keep_by_loc: dict[str, set] = {}
    for loc, name in keep:
        loc = normalize_location(loc)
        keep_by_loc.setdefault(loc, set())
        if name in location_map(mesh, loc):
            keep_by_loc[loc].add(name)
        elif not ignore_missing:
            require_key(mesh, loc, name)

    drop_by_loc: dict[str, set] = {loc: set() for loc in LOCATIONS}
    for loc, name in drop:
        loc = normalize_location(loc)
        if name in location_map(mesh, loc):
            drop_by_loc[loc].add(name)
        elif not ignore_missing:
            require_key(mesh, loc, name)

    rename_by_loc: dict[str, dict] = {loc: {} for loc in LOCATIONS}
    targets: dict[str, set] = {loc: set() for loc in LOCATIONS}
    for loc, old, new in rename:
        loc = normalize_location(loc)
        if old not in location_map(mesh, loc):
            if ignore_missing:
                continue
            require_key(mesh, loc, old)
        if not new:
            raise ValueError(
                f"meshio++: data_manage: cannot rename {loc} '{old}' to an empty name"
            )
        if old in rename_by_loc[loc]:
            raise ValueError(f"meshio++: data_manage: {loc} '{old}' is renamed twice")
        if new in targets[loc]:
            raise ValueError(
                f"meshio++: data_manage: two renames both target {loc} '{new}'"
            )
        targets[loc].add(new)
        rename_by_loc[loc][old] = new

    # A rename whose target already exists would silently clobber it, unless
    # that target is itself renamed away or dropped in the same pass.
    for loc, renames in rename_by_loc.items():
        for old, new in renames.items():
            if new == old:
                continue
            exists = new in location_map(mesh, loc)
            renamed_away = new in renames
            dropped = new in drop_by_loc[loc] or (
                loc in keep_by_loc and new not in keep_by_loc[loc]
            )
            if exists and not renamed_away and not dropped:
                raise ValueError(
                    f"meshio++: data_manage: cannot rename {loc} '{old}' to "
                    f"'{new}' (that name already exists)"
                )

    out = _clone(mesh)
    dropped_report = []
    renamed_report = []
    for loc in LOCATIONS:
        source = location_map(mesh, loc)
        result = {}
        for name in sorted(source.keys()):
            if loc in keep_by_loc and name not in keep_by_loc[loc]:
                dropped_report.append(f"{loc}:{name}")
                continue
            if name in drop_by_loc[loc]:
                dropped_report.append(f"{loc}:{name}")
                continue
            target = rename_by_loc[loc].get(name, name)
            if target != name:
                renamed_report.append((f"{loc}:{name}", f"{loc}:{target}"))
            result[target] = copy.deepcopy(source[name])
        setattr(out, loc, result)

    return {
        "mesh": out,
        "dropped": sorted(dropped_report),
        "renamed": renamed_report,
    }


def data_manage(mesh, keep=None, drop=None, rename=None, ignore_missing=False) -> dict:
    """Rewrite which data arrays ``mesh`` carries.

    Args:
        mesh: the mesh to rewrite (unmodified).
        keep: iterable of ``(location, name)`` to retain. Only affects the
            locations it mentions.
        drop: iterable of ``(location, name)`` to remove, applied after ``keep``.
        rename: iterable of ``(location, old, new)``, applied last.
        ignore_missing: skip keys that do not exist instead of raising.

    Returns:
        ``{"mesh": Mesh, "dropped": [str], "renamed": [(str, str)]}`` where the
        report entries are qualified as ``"point_data:T"``.

    Raises:
        ValueError: an unknown key (unless ``ignore_missing``), or a rename
            collision.
    """
    keep = list(keep or [])
    drop = list(drop or [])
    rename = [tuple(r) for r in (rename or [])]

    try:
        from . import _core
    except Exception:
        return _manage_py(mesh, keep, drop, rename, ignore_missing)

    # A ValueError is the user's own error (unknown key, rename collision) and
    # must propagate; anything else means the core could not handle this mesh,
    # so fall back to the reference implementation.
    try:
        raw = _core.data_manage(mesh, keep, drop, rename, ignore_missing)
    except ValueError:
        raise
    except Exception:
        return _manage_py(mesh, keep, drop, rename, ignore_missing)

    out = raw["mesh"]
    _carry_sets(mesh, out)
    return {
        "mesh": out,
        "dropped": list(raw["dropped"]),
        "renamed": [tuple(t) for t in raw["renamed"]],
    }


def _carry_sets(source, target) -> None:
    """Copy ``point_sets``/``cell_sets`` across, which the C++ core never sees.

    They index geometry, and the data operations never modify geometry, so they
    survive verbatim.
    """
    for attr in ("point_sets", "cell_sets"):
        value = getattr(source, attr, None)
        if value:
            setattr(target, attr, copy.deepcopy(value))


def data_drop(mesh, location: str, names, ignore_missing: bool = False):
    """Drop the named arrays at one location.

    Args:
        mesh: the mesh to rewrite (unmodified).
        location: ``"point"``, ``"cell"`` or ``"field"``.
        names: iterable of array names to remove.
        ignore_missing: skip names that do not exist instead of raising.

    Returns:
        The rewritten mesh.
    """
    loc = normalize_location(location)
    return data_manage(
        mesh, drop=[(loc, n) for n in names], ignore_missing=ignore_missing
    )["mesh"]


def data_keep(mesh, location: str, names, ignore_missing: bool = False):
    """Keep only the named arrays at one location, dropping the rest there.

    The other two locations are left untouched.

    Args:
        mesh: the mesh to rewrite (unmodified).
        location: ``"point"``, ``"cell"`` or ``"field"``.
        names: iterable of array names to retain.
        ignore_missing: skip names that do not exist instead of raising.

    Returns:
        The rewritten mesh.
    """
    loc = normalize_location(location)
    wanted = set()
    for n in names:
        if n in location_map(mesh, loc):
            wanted.add(n)
        elif not ignore_missing:
            require_key(mesh, loc, n)
    # Expressed as the complement, so an empty `names` means "drop the lot"
    # without needing a "whitelist is active" sentinel.
    unwanted = [n for n in location_map(mesh, loc) if n not in wanted]
    return data_manage(mesh, drop=[(loc, n) for n in unwanted])["mesh"]


def data_rename(mesh, location: str, old: str, new: str):
    """Rename one array, preserving its values, dtype and shape.

    Args:
        mesh: the mesh to rewrite (unmodified).
        location: ``"point"``, ``"cell"`` or ``"field"``.
        old: the existing name.
        new: the new name.

    Returns:
        The rewritten mesh.
    """
    return data_manage(mesh, rename=[(normalize_location(location), old, new)])["mesh"]
