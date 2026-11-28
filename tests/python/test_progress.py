"""Tests for cooperative cancellation and progress reporting (roadmap 9.16)."""

from __future__ import annotations

import pytest

import meshioplusplus as mio


def test_progress_context_reports_octree_stages():
    box = mio.box((2, 2, 2), cell_type="tetra")
    skin = mio.extract_skin(box)
    events = []

    def on_progress(fraction: float, stage: str) -> bool:
        events.append((fraction, stage))
        return True

    with mio.progress_context(callback=on_progress):
        res = mio.octree_mesh(skin, min_size=0.5, max_size=0.5)

    assert len(res.cells) > 0
    assert len(events) >= 3
    assert events[0][1] == "octree_mesh:prepare"
    assert events[-1][1] == "octree_mesh:done"


def test_progress_callback_returning_false_cancels():
    box = mio.box((2, 2, 2), cell_type="tetra")
    skin = mio.extract_skin(box)
    seen = []

    def cancel_on_second(fraction: float, stage: str) -> bool:
        seen.append(stage)
        return len(seen) < 2

    with mio.progress_context(callback=cancel_on_second):
        with pytest.raises(mio.CancelledError):
            mio.octree_mesh(skin, min_size=0.5, max_size=0.5)

    assert len(seen) == 2


def test_cancel_token_cancels_adapt_and_optimize():
    box = mio.box((2, 2, 2), cell_type="tetra")
    token = mio.CancelToken()
    assert not token.is_cancelled
    token.cancel()
    assert token.is_cancelled

    with mio.progress_context(cancel_token=token):
        with pytest.raises(mio.CancelledError):
            mio.optimize_volume(box, max_iterations=3)
        with pytest.raises(mio.CancelledError):
            mio.adapt_volume(box, hmax=0.4, max_iterations=3)

    token.reset()
    assert not token.is_cancelled


def test_progress_callback_exception_propagates():
    box = mio.box((2, 2, 2), cell_type="tetra")
    skin = mio.extract_skin(box)

    def boom(fraction: float, stage: str) -> bool:
        raise RuntimeError("boom from callback")

    with mio.progress_context(callback=boom):
        with pytest.raises(RuntimeError, match="boom from callback"):
            mio.octree_mesh(skin, min_size=0.5, max_size=0.5)
