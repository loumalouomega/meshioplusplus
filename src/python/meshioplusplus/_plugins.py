"""Discovery of the add-on wheels (``meshioplusplus-occt`` and its kin).

An add-on registers its provider with the shared core when its extension is
imported. The base package finds add-ons through the ``meshioplusplus.plugins``
entry-point group, once and lazily, from the capability queries and the paths
that need a provider. An add-on must come from the same release as the base
(its wheel pins it exactly); a mismatch is refused by name rather than loaded.
``MESHIOPLUSPLUS_NO_PLUGINS=1`` skips discovery, so a test run can prove the
base wheel alone.
"""

import os
import threading
import warnings
from importlib import import_module

GROUP = "meshioplusplus.plugins"

_lock = threading.Lock()
_state = {"done": False, "loaded": {}, "errors": {}}


def _entry_points():
    try:
        from importlib.metadata import entry_points
    except ImportError:  # pragma: no cover
        return []
    try:
        eps = entry_points()
        if hasattr(eps, "select"):
            return list(eps.select(group=GROUP))
        return list(eps.get(GROUP, []))  # Python < 3.10
    except Exception:  # a broken metadata directory must not break the import
        return []


def _dist_version(entry_point):
    """The version of the distribution that supplies ``entry_point``, or None.

    Read from the package metadata, without importing the add-on: importing it
    registers its provider with the core, so the check must come first.
    """
    dist = getattr(entry_point, "dist", None)
    if dist is not None:
        return dist.version
    try:  # Python < 3.10: entry points carry no distribution
        from importlib.metadata import distributions
    except ImportError:  # pragma: no cover
        return None
    for candidate in distributions():
        for ep in candidate.entry_points:
            if (
                ep.group == GROUP
                and ep.name == entry_point.name
                and ep.value == entry_point.value
            ):
                return candidate.version
    return None


def _base_version():
    from .__about__ import __version__

    return __version__


def load_plugins(force=False):
    """Import every installed add-on once; return ``{entry name: module}``."""
    with _lock:
        if _state["done"] and not force:
            return dict(_state["loaded"])
        _state["done"] = True
        _state["loaded"], _state["errors"] = {}, {}
        if os.environ.get("MESHIOPLUSPLUS_NO_PLUGINS", "") not in ("", "0"):
            return {}
        base = _base_version()
        for ep in _entry_points():
            try:
                # An add-on of another release is refused before it is imported:
                # importing it registers its provider with the core.
                dist_version = _dist_version(ep)
                if dist_version is not None and dist_version != base:
                    raise ImportError(
                        f"add-on {ep.name!r} is version {dist_version}, but meshioplusplus is "
                        f"{base}; install the matching release (pip install -U "
                        f"meshioplusplus[{ep.name}])"
                    )
                module = import_module(ep.value.split(":")[0].strip())
                version = getattr(module, "__version__", None)
                if version != base:
                    raise ImportError(
                        f"add-on {ep.name!r} is version {version}, but meshioplusplus is "
                        f"{base}; install the matching release (pip install -U "
                        f"meshioplusplus[{ep.name}])"
                    )
                _state["loaded"][ep.name] = module
            except Exception as exc:
                _state["errors"][ep.name] = str(exc)
                warnings.warn(f"meshioplusplus: add-on {ep.name!r} not loaded: {exc}")
        return dict(_state["loaded"])


def plugin_errors():
    """Why each add-on that failed to load was refused."""
    load_plugins()
    return dict(_state["errors"])
