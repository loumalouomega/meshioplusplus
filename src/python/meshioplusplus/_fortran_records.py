"""Split a Fortran sequential unformatted file into its records.

The Python twin of ``detail/fortran_records.hpp``: each record is framed by its
byte length before and after it, in a 4- or 8-byte marker stored in the
writing machine's byte order. Abaqus ``.fil``, Nastran OP2 and the OpenRadioss
time-history file share it; the last one also writes through it
(:func:`fortran_record`).
"""

import struct

from ._exceptions import ReadError, WriteError

__all__ = ["sniff_fortran_records", "fortran_records", "fortran_record"]


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


def fortran_records(data, layout, what, lenient=False):
    """Every record's payload span ``(offset, size)``, in file order.

    With ``lenient`` the first record that runs past the end of the buffer, or
    whose trailing marker does not match its leading one, ends the list instead
    of raising (a solver killed mid-write leaves such a tail).
    """
    width, order = layout
    fmt = order + ("i" if width == 4 else "q")
    out = []
    pos = 0
    size = len(data)
    while pos < size:
        if size - pos < width:
            if lenient:
                break
            raise ReadError(
                f"{what}: Fortran record at offset {pos} runs past the end of the file"
            )
        n = struct.unpack(fmt, data[pos : pos + width])[0]
        if n < 0 or n > size - pos - width:
            if lenient:
                break
            raise ReadError(
                f"{what}: Fortran record at offset {pos} runs past the end of the file"
            )
        tail = data[pos + width + n : pos + 2 * width + n]
        if len(tail) != width or struct.unpack(fmt, tail)[0] != n:
            if lenient:
                break
            raise ReadError(
                f"{what}: Fortran record at offset {pos} has mismatched length markers"
            )
        out.append((pos + width, n))
        pos += 2 * width + n
    return out


def fortran_record(payload, layout=(4, ">")):
    """``payload`` framed as one record: its length, the bytes, the length again.

    ``layout`` is ``(marker bytes, byte order)`` as :func:`sniff_fortran_records`
    returns it. The inverse of :func:`fortran_records`; the C++ twin is
    ``fortran_append_record``.
    """
    width, order = layout
    if width == 4 and len(payload) > 0x7FFFFFFF:
        raise WriteError(
            f"Fortran records: a record of {len(payload)} bytes does not fit a "
            "4-byte length marker"
        )
    marker = struct.pack(order + ("i" if width == 4 else "q"), len(payload))
    return marker + bytes(payload) + marker
