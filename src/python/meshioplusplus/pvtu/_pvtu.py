"""VTK XML parallel index ``.pvtu``, the pure-Python reference.

A ``<PUnstructuredGrid>`` index holds the *declarations* (``PPoints``, ``PPointData``,
``PCellData``) and one ``<Piece Source=>`` per part; every piece is a
standalone ``.vtu`` file. The machinery is shared with the other index
formats and lives in ``meshioplusplus._pvtk_index``.
"""

from __future__ import annotations

from .. import _pvtk_index as _ix

KIND = "pvtu"


def read(filename, piece=None, ghosts="keep"):
    return _ix.read_index(filename, KIND, piece=piece, ghosts=ghosts)


def write(
    filename,
    mesh,
    binary=True,
    compression="zlib",
    header_type=None,
    part_key=_ix.PARTITION_PART,
):
    _ix.write_mesh(
        filename,
        mesh,
        KIND,
        part_key=part_key,
        binary=binary,
        compression=compression,
        header_type=header_type,
    )


def write_pieces(filename, pieces, binary=True, compression="zlib", header_type=None):
    _ix.write_pieces(
        filename,
        pieces,
        KIND,
        binary=binary,
        compression=compression,
        header_type=header_type,
    )
