"""A ZIP writer that matches ``detail/zip_archive.hpp`` byte for byte.

Members are deflated at level 9 (raw deflate, the stream ``gzip_compress`` makes)
when that is smaller, else stored; every member is dated 1980-01-01 00:00 with no
extra fields, so the archive depends on the contents alone. Reading goes through
``zipfile``. No Zip64.
"""

from __future__ import annotations

import struct
import zlib

from ._exceptions import WriteError


class ZipWriter:
    def __init__(self):
        self._out = bytearray()
        self._directory = bytearray()
        self._count = 0

    def add(self, name, data, deflate=True):
        if isinstance(data, str):
            data = data.encode("utf-8")
        raw_name = name.encode("utf-8")
        body, method = data, 0
        if deflate and data:
            c = zlib.compressobj(9, zlib.DEFLATED, -15, 8, zlib.Z_DEFAULT_STRATEGY)
            packed = c.compress(data) + c.flush()
            if len(packed) < len(data):
                body, method = packed, 8
        crc = zlib.crc32(data) & 0xFFFFFFFF
        offset = len(self._out)
        if (
            offset + len(body) > 0xFFFFFFFF
            or len(data) > 0xFFFFFFFF
            or self._count == 0xFFFF
        ):
            raise WriteError(
                "meshio++: zip: the archive would need Zip64, which is not written"
            )
        common = struct.pack(
            "<HHHHHIIIHH",
            20,
            0,
            method,
            0,
            0x21,
            crc,
            len(body),
            len(data),
            len(raw_name),
            0,
        )
        self._out += struct.pack("<I", 0x04034B50) + common + raw_name + body
        self._directory += (
            struct.pack("<IH", 0x02014B50, 20)
            + common
            + struct.pack("<HHHII", 0, 0, 0, 0, offset)
            + raw_name
        )
        self._count += 1

    def finish(self):
        out = bytes(self._out) + bytes(self._directory)
        if len(out) > 0xFFFFFFFF:
            raise WriteError(
                "meshio++: zip: the archive would need Zip64, which is not written"
            )
        return out + struct.pack(
            "<IHHHHIIH",
            0x06054B50,
            0,
            0,
            self._count,
            self._count,
            len(self._directory),
            len(self._out),
            0,
        )
