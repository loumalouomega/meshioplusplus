"""A PNG encoder for RGBA8 (or RGB8) images, for the texture options (roadmap
2.4.8): rows top to bottom, filter 0, one zlib stream."""

from __future__ import annotations

import os
import struct
import zlib

import numpy as np

SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _chunk(kind, data):
    body = kind + data
    return (
        struct.pack(">I", len(data))
        + body
        + struct.pack(">I", zlib.crc32(body) & 0xFFFFFFFF)
    )


def encode(image) -> bytes:
    """``image`` ``(height, width, 3 or 4)`` uint8 as PNG bytes."""
    a = np.asarray(image)
    if a.ndim != 3 or a.shape[2] not in (3, 4) or a.shape[0] < 1 or a.shape[1] < 1:
        raise ValueError(
            "meshio++: a texture image must be (height, width, 3 or 4), got "
            f"{a.shape}"
        )
    a = np.ascontiguousarray(np.clip(a, 0, 255).astype(np.uint8))
    h, w, c = a.shape
    raw = b"".join(b"\0" + a[r].tobytes() for r in range(h))
    header = struct.pack(">IIBBBBB", w, h, 8, 6 if c == 4 else 2, 0, 0, 0)
    return (
        SIGNATURE
        + _chunk(b"IHDR", header)
        + _chunk(b"IDAT", zlib.compress(raw, 9))
        + _chunk(b"IEND", b"")
    )


def texture_bytes(texture) -> bytes:
    """A texture option as PNG bytes: a ``.png`` path, PNG bytes, or an image array."""
    if isinstance(texture, (bytes, bytearray)):
        data = bytes(texture)
    elif isinstance(texture, (str, os.PathLike)):
        with open(texture, "rb") as f:
            data = f.read()
    else:
        return encode(texture)
    if not data.startswith(SIGNATURE):
        raise ValueError("meshio++: a texture file must be a PNG image")
    return data
