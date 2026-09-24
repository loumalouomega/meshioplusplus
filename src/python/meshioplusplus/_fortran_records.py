"""Split a Fortran sequential unformatted file into its records.

The Python twin of ``detail/fortran_records.hpp``: each record is framed by its
byte length before and after it, in a 4- or 8-byte marker stored in the
writing machine's byte order. Abaqus ``.fil`` and Nastran OP2 share it.
"""

import struct

from ._exceptions import ReadError

__all__ = ["sniff_fortran_records", "fortran_records"]


def sniff_fortran_records(data):
    """(marker bytes, byte order) framing the first record, or ``None``.

    Tries 4-byte then 8-byte markers, each little- then big-endian; a candidate
    is accepted when its length is positive, fits the buffer and the same value
    follows the payload.
    """
    for width in (4, 8):
        for order in ("<", ">"):
            if len(data) < width:
                continue
            fmt = order + ("i" if width == 4 else "q")
            n = struct.unpack(fmt, data[:width])[0]
            if n <= 0 or n > len(data) - width:
                continue
            tail = data[width + n : 2 * width + n]
            if len(tail) == width and struct.unpack(fmt, tail)[0] == n:
                return width, order
    return None


def fortran_records(data, layout, what):
    """Every record's payload span ``(offset, size)``, in file order."""
    width, order = layout
    fmt = order + ("i" if width == 4 else "q")
    out = []
    pos = 0
    size = len(data)
    while pos < size:
        if size - pos < width:
            raise ReadError(
                f"{what}: Fortran record at offset {pos} runs past the end of the file"
            )
        n = struct.unpack(fmt, data[pos : pos + width])[0]
        if n < 0 or n > size - pos - width:
            raise ReadError(
                f"{what}: Fortran record at offset {pos} runs past the end of the file"
            )
        tail = data[pos + width + n : pos + 2 * width + n]
        if len(tail) != width or struct.unpack(fmt, tail)[0] != n:
            raise ReadError(
                f"{what}: Fortran record at offset {pos} has mismatched length markers"
            )
        out.append((pos + width, n))
        pos += 2 * width + n
    return out
