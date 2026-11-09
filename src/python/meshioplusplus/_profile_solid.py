"""Profile solids through a CAD provider (roadmap 7.7.4): a closed
``fit_curves`` sketch extruded or revolved, then filleted, chamfered or
shelled. C++ core and an OpenCASCADE provider only (see
``doc/profile_solid.md``)."""

from __future__ import annotations

__all__ = ["profile_solid", "profile_solid_providers", "profile_features_from_text"]

_KEYS = (
    "kind",
    "distance",
    "symmetric",
    "angle",
    "axis_origin",
    "axis_direction",
    "radius",
)


def profile_solid(
    sketch,
    features,
    step_path: str | None = None,
    *,
    provider: str = "",
    deflection: float = 0.0,
    scan=None,
    return_report: bool = False,
):
    """Build a solid from a ``fit_curves`` sketch by a list of features
    applied in order: the first ``{"kind": "extrude", "distance": d}`` (along
    the sketch's normal; ``0`` runs to ``scan``'s extent; ``"symmetric"``
    half each way) or ``{"kind": "revolve", "axis_origin": p,
    "axis_direction": d, "angle": 360}``, then ``{"kind": "fillet",
    "radius": r}``, ``{"kind": "chamfer", "distance": d}`` (every edge) or
    ``{"kind": "shell", "distance": wall}`` (the end face removed).

    :returns: the solid's tessellation; with ``return_report`` ``(mesh,
        report)``: ``volume``, ``area``, ``num_faces``, ``valid``,
        ``provider`` and ``warnings``. ``step_path`` receives it as STEP.
    :raises ValueError: on a mesh without a sketch, loops that do not close or
        a bad feature.
    :raises NotImplementedError: without the C++ core or a profile-solid
        provider (OpenCASCADE: a build with it, or ``meshioplusplus-occt``).
    """
    from ._hull import _core_op

    if isinstance(features, str):
        features = profile_features_from_text(features)
    items = []
    for f in features:
        unknown = set(f) - set(_KEYS)
        if unknown:
            raise ValueError(
                f"meshio++: profile_solid: a feature has no key {sorted(unknown)[0]!r}"
            )
        item = dict(f)
        for key in ("axis_origin", "axis_direction"):
            if key in item:
                item[key] = [float(v) for v in item[key]]
        items.append(item)
    res = _core_op("profile_solid")(
        sketch, items, str(step_path or ""), str(provider), float(deflection), scan
    )
    out = res.pop("mesh")
    return (out, res) if return_report else out


def profile_solid_providers():
    """The names of the registered profile-solid providers (``"occt"`` in a
    build with OpenCASCADE); empty without one."""
    try:
        from . import _core
    except ImportError:
        return []
    names = getattr(_core, "profile_solid_provider_names", None)
    return list(names()) if names else []


def profile_features_from_text(text):
    """Features from their text form, ``;``-separated, as the CLIs and the
    pipeline take them: ``extrude:D`` (``extrude:D:symmetric``),
    ``revolve:OX,OY,OZ:DX,DY,DZ[:ANGLE]``, ``fillet:R``, ``chamfer:D``,
    ``shell:W``."""

    def bad(item, why):
        return ValueError(f"meshio++: profile_solid: {why}, not {item!r}")

    def vec3(item, part):
        c = part.split(",")
        if len(c) != 3:
            raise bad(item, f"{part!r} is not x,y,z")
        return [float(v) for v in c]

    out = []
    for item in filter(None, text.split(";")):
        parts = item.split(":")
        kind = parts[0]
        if kind == "extrude":
            if len(parts) not in (2, 3) or (
                len(parts) == 3 and parts[2] != "symmetric"
            ):
                raise bad(item, "an extrusion is extrude:D or extrude:D:symmetric")
            out.append(
                {
                    "kind": kind,
                    "distance": float(parts[1]),
                    "symmetric": len(parts) == 3,
                }
            )
        elif kind == "revolve":
            if len(parts) not in (3, 4):
                raise bad(item, "a revolution is revolve:OX,OY,OZ:DX,DY,DZ[:ANGLE]")
            f = {
                "kind": kind,
                "axis_origin": vec3(item, parts[1]),
                "axis_direction": vec3(item, parts[2]),
            }
            if len(parts) == 4:
                f["angle"] = float(parts[3])
            out.append(f)
        elif kind == "fillet":
            if len(parts) != 2:
                raise bad(item, "a fillet is fillet:R")
            out.append({"kind": kind, "radius": float(parts[1])})
        elif kind in ("chamfer", "shell"):
            if len(parts) != 2:
                raise bad(item, f"a {kind} takes one distance")
            out.append({"kind": kind, "distance": float(parts[1])})
        else:
            raise bad(item, "a feature is extrude, revolve, fillet, chamfer or shell")
    return out
