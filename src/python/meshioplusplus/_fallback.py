"""Deciding whether a declined C++ fast path may fall back to its Python twin.

Every format package wraps ``_core.<fmt>_read`` / ``_write`` in a shim that
retries the pure-Python reference implementation when the C++ core cannot
handle a file. That fallback used to be a bare ``except Exception: pass``,
which hid *why* the fast path declined and swallowed bugs (a parser tripping
over ``std::stoll``) along with legitimate declines (a construct the core
deliberately does not handle). :func:`core_declined` is the one place that
decides, and says so.

======================================  =========================  ==========
exception out of ``_core``              meaning                    default
======================================  =========================  ==========
``ReadError`` / ``WriteError``          a recognised decline       DEBUG, fall back
``TypeError`` / ``MemoryError`` /       never a decline            propagate
``RecursionError``
anything else                           the fast path *broke*      WARNING, fall back
======================================  =========================  ==========

The last row is deliberate: the C++ readers throw ``ReadError`` for nearly every
refusal, but a bare ``std::stoll``/``.at()`` on a malformed token surfaces as
``ValueError``/``IndexError``, and an ambiguous ``.msh`` relies on the Python
twin turning that into a clean ``ReadError`` so the next candidate is tried.
Propagating it would trade a silent fallback for a hard failure; a warning
makes it visible without that regression. The same goes for an
``AttributeError``: a shim naming a ``_core`` entry point the compiled extension
does not have (a stale editable build, or a format compiled out) still works
through the Python twin, but says so instead of silently running slower.

Set ``MESHIOPLUSPLUS_STRICT_CORE=1`` (also ``on``/``true``/``yes``) to turn every
fall-back into a re-raise. That is how a contributor finds out that a format
believed to be native has been quietly running its Python twin.
"""

from __future__ import annotations

import logging
import os

from ._exceptions import ReadError, WriteError

_LOG = logging.getLogger("meshioplusplus")

_ENV_VAR = "MESHIOPLUSPLUS_STRICT_CORE"
_TRUTHY = ("1", "on", "true", "yes")

# `RecursionError` subclasses `RuntimeError`, so it has to be named here to
# stay out of the "fast path broke, warn and fall back" bucket.
_NEVER_A_DECLINE = (TypeError, MemoryError, RecursionError)
_RECOGNISED_DECLINE = (ReadError, WriteError)

_strict: bool | None = None


def strict_core() -> bool:
    """Whether ``MESHIOPLUSPLUS_STRICT_CORE`` is on (read once, then cached)."""
    global _strict
    if _strict is None:
        _strict = os.environ.get(_ENV_VAR, "").strip().lower() in _TRUTHY
    return _strict


def set_strict_core(enabled: bool | None) -> None:
    """Override the environment switch; ``None`` re-reads it on next use."""
    global _strict
    _strict = enabled


def core_declined(exc: BaseException, fmt: str, direction: str, target) -> bool:
    """Log why the C++ fast path failed and say whether to fall back.

    Call it from the shim's ``except Exception as exc`` block::

        except Exception as exc:
            if not core_declined(exc, "ansys", "read", filename):
                raise

    :param exc: the exception the ``_core`` call raised.
    :param fmt: the format name, for the log line.
    :param direction: ``"read"`` or ``"write"``.
    :param target: the path being read or written, for the log line.
    :returns: ``True`` to fall back to the Python twin, ``False`` to re-raise.
    """
    if isinstance(exc, _NEVER_A_DECLINE):
        return False

    if strict_core():
        _LOG.warning(
            "meshio++: %s %s '%s': C++ core refused (%s: %s); %s=1, not falling back",
            fmt,
            direction,
            target,
            type(exc).__name__,
            exc,
            _ENV_VAR,
        )
        return False

    if isinstance(exc, _RECOGNISED_DECLINE):
        _LOG.debug(
            "meshio++: %s %s '%s': C++ core declined (%s: %s); using the Python twin",
            fmt,
            direction,
            target,
            type(exc).__name__,
            exc,
        )
    else:
        _LOG.warning(
            "meshio++: %s %s '%s': C++ core failed unexpectedly (%s: %s); "
            "using the Python twin. Set %s=1 to raise instead.",
            fmt,
            direction,
            target,
            type(exc).__name__,
            exc,
            _ENV_VAR,
        )
    return True
