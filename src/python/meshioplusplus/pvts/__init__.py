"""VTK XML parallel index ``.pvts``: a whole structured grid placed from its pieces. Python only, read-only."""

from .. import _pvtk_structured as _ps
from .._helpers import register_format


def read(filename, piece=None):
    """Read a ``.pvts`` index: every piece placed by its ``Extent`` in ``WholeExtent``
    (shared planes become one layer of points), one ``piece_<i>`` cell region each, or
    one piece with ``piece=``."""
    return _ps.read(filename, "pvts", piece=piece)


register_format("pvts", [".pvts"], read, {})

__all__ = ["read"]
