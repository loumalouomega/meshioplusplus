"""VTK XML parallel index ``.pvti``: a whole structured grid placed from its pieces. Python only, read-only."""

from .. import _pvtk_structured as _ps
from .._helpers import register_format


def read(filename, piece=None):
    """Read a ``.pvti`` index: every piece placed by its ``Extent`` in ``WholeExtent``
    (shared planes become one layer of points), one ``piece_<i>`` cell region each, or
    one piece with ``piece=``."""
    return _ps.read(filename, "pvti", piece=piece)


register_format("pvti", [".pvti"], read, {})

__all__ = ["read"]
