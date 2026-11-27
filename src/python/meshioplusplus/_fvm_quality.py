"""Finite-volume (CFD) interior face-pair quality diagnostics (roadmap 9.11).

C++-core only. Computes non-orthogonality, FVM face skewness, adjacent volume
growth ratio, and face pyramid-height checks across interior face pairs.
"""

from __future__ import annotations

__all__ = [
    "compute_fvm_quality",
    "attach_fvm_quality",
    "FVM_METRIC_NAMES",
]

FVM_METRIC_NAMES = (
    "fvm:non_orthogonality",
    "fvm:face_skewness",
    "fvm:volume_ratio",
    "fvm:min_pyramid_height",
)


def compute_fvm_quality(mesh):
    """Compute finite-volume face-pair quality diagnostics on a 3-D volume mesh."""
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "compute_fvm_quality"):
        raise NotImplementedError(
            "meshio++: compute_fvm_quality: requires the compiled meshioplusplus._core extension"
        )
    return _core.compute_fvm_quality(mesh)


def attach_fvm_quality(mesh):
    """Return a copy of ``mesh`` with ``fvm:*`` cell_data arrays attached."""
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None or not hasattr(_core, "attach_fvm_quality"):
        raise NotImplementedError(
            "meshio++: attach_fvm_quality: requires the compiled meshioplusplus._core extension"
        )
    return _core.attach_fvm_quality(mesh)
