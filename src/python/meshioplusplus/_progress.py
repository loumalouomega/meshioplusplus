"""Cooperative cancellation and progress reporting helpers (roadmap 9.16)."""

from __future__ import annotations

from contextlib import contextmanager
from typing import Callable, Iterator, Optional

ProgressCallback = Callable[[float, str], Optional[bool]]


class CancelToken:
    """Thread-safe cancellation token for long-running C++ operations."""

    def __init__(self) -> None:
        try:
            from . import _core

            self._token = _core.CancelToken()
        except ImportError:
            self._token = None
            self._flag = False

    def cancel(self) -> None:
        if self._token is not None:
            self._token.cancel()
        else:
            self._flag = True

    def reset(self) -> None:
        if self._token is not None:
            self._token.reset()
        else:
            self._flag = False

    @property
    def is_cancelled(self) -> bool:
        if self._token is not None:
            return bool(self._token.is_cancelled)
        return bool(self._flag)


@contextmanager
def progress_context(
    callback: Optional[ProgressCallback] = None,
    cancel_token: Optional[CancelToken] = None,
) -> Iterator[None]:
    """Install a thread-local progress callback and/or cancel token.

    ``callback(fraction, stage)`` receives a float in ``[0.0, 1.0]`` (or ``< 0``
    when indeterminate) and a stage label string; returning ``False`` cancels
    the active operation with :class:`CancelledError`.
    """
    if callback is None and cancel_token is None:
        yield
        return
    try:
        from . import _core
    except ImportError:
        yield
        return
    raw_token = cancel_token._token if cancel_token is not None else None
    prev = _core._push_progress_context(callback, raw_token)
    try:
        yield
    finally:
        _core._pop_progress_context(prev)
