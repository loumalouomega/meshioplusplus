"""MMG remesher for meshio++: importing it registers the provider."""

from . import _mmg

__version__ = _mmg.__version__
__mmg_version__ = _mmg.__mmg_version__

__all__ = ["__version__", "__mmg_version__"]
