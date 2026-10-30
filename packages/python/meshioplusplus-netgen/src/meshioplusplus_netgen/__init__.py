"""Netgen tetrahedral mesher for meshio++: importing it registers the provider."""

from . import _netgen

__version__ = _netgen.__version__
__netgen_version__ = _netgen.__netgen_version__

__all__ = ["__version__", "__netgen_version__"]
