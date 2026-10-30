"""OpenCASCADE add-on for meshio++: importing it registers the CAD kernel."""

from . import _occt

__version__ = _occt.__version__
__occt_version__ = _occt.__occt_version__

__all__ = ["__version__", "__occt_version__"]
