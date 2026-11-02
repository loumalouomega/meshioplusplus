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


# --------------------------------------------------------------------------- #
# Third-party formats (roadmap 7.1.3)                                          #
# --------------------------------------------------------------------------- #
#: The entry-point group a separate package declares a format in. Its value is
#: either a callable taking no arguments, which registers the package's formats
#: (``meshioplusplus.register_format``), or a module that registers them when
#: imported. Unlike the add-ons, a format package is not pinned to a release:
#: ``register_format`` is the stable contract.
FORMATS_GROUP = "meshioplusplus.formats"

_format_lock = threading.RLock()
_format_state = {"done": False, "loading": False, "loaded": {}, "errors": {}}


def _format_entry_points():
    try:
        from importlib.metadata import entry_points
    except ImportError:  # pragma: no cover
        return []
    try:
        eps = entry_points()
        if hasattr(eps, "select"):
            return list(eps.select(group=FORMATS_GROUP))
        return list(eps.get(FORMATS_GROUP, []))  # pragma: no cover
    except Exception:  # a broken metadata directory must not break a read
        return []


def load_format_plugins(force=False):
    """Register every installed third-party format once; return
    ``{entry name: loaded object}``.

    Called on the first use of the format registry (a read, a write, a
    registry listing), so installing a package that declares a
    ``meshioplusplus.formats`` entry point is all it takes. A plugin that fails
    to load is skipped with a warning naming it, and ``format_plugin_errors``
    says why; ``MESHIOPLUSPLUS_NO_PLUGINS=1`` skips discovery, as for add-ons.
    """
    with _format_lock:
        # `loading`: a plugin that reads the registry while it registers comes
        # back here on the same thread (the lock is reentrant); other threads
        # wait for the lock and then see `done`.
        if _format_state["loading"] or (_format_state["done"] and not force):
            return dict(_format_state["loaded"])
        _format_state["loaded"], _format_state["errors"] = {}, {}
        if os.environ.get("MESHIOPLUSPLUS_NO_PLUGINS", "") not in ("", "0"):
            _format_state["done"] = True
            return {}
        _format_state["loading"] = True
        try:
            _load_format_entry_points()
        finally:
            _format_state["loading"] = False
            _format_state["done"] = True
        return dict(_format_state["loaded"])


def _load_format_entry_points():
    for ep in sorted(_format_entry_points(), key=lambda e: e.name):
        try:
            obj = ep.load()
            if callable(obj):
                obj()
            _format_state["loaded"][ep.name] = obj
        except Exception as exc:
            _format_state["errors"][ep.name] = f"{type(exc).__name__}: {exc}"
            warnings.warn(
                f"meshioplusplus: format plugin {ep.name!r} ({ep.value}) not "
                f"loaded: {exc}",
                RuntimeWarning,
                stacklevel=2,
            )


def format_plugin_errors():
    """Why each format plugin that failed to load was skipped."""
    load_format_plugins()
    return dict(_format_state["errors"])
