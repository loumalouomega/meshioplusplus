"""
LZF stream codec for the PCD ``binary_compressed`` payload.

The stream format is the one liblzf defines (<http://oldhome.schmorp.de/marc/liblzf.html>):
a control byte below 32 introduces a literal run of ``ctrl + 1`` bytes; anything
else is a back reference whose length lives in the top three bits (7 meaning "add
the next byte") and whose 13-bit distance is split over the low five bits and one
following byte. This is an independent implementation of that stream format -- the
decoder follows it exactly, and the encoder emits a valid (not byte-identical to
liblzf) stream: it uses a single-slot 3-byte-hash table, like liblzf, but only the
decoder's output is normative.
"""

from .._exceptions import ReadError

_MAX_LITERAL = 32
_MAX_REF = 264  # 2 + 7 + 255
_MAX_OFFSET = 1 << 13
_HASH_BITS = 14


def decompress(data, out_size):
    """Decode an LZF stream into exactly ``out_size`` bytes."""
    src = bytes(data)
    n = len(src)
    out = bytearray(out_size)
    ip = 0
    op = 0
    while ip < n:
        ctrl = src[ip]
        ip += 1
        if ctrl < 32:
            length = ctrl + 1
            if ip + length > n:
                raise ReadError("PCD: truncated LZF literal run")
            if op + length > out_size:
                raise ReadError("PCD: LZF stream overruns the declared size")
            out[op : op + length] = src[ip : ip + length]
            ip += length
            op += length
        else:
            length = ctrl >> 5
            if length == 7:
                if ip >= n:
                    raise ReadError("PCD: truncated LZF back reference")
                length += src[ip]
                ip += 1
            if ip >= n:
                raise ReadError("PCD: truncated LZF back reference")
            ref = op - ((ctrl & 0x1F) << 8) - src[ip] - 1
            ip += 1
            length += 2
            if ref < 0:
                raise ReadError("PCD: LZF back reference before the start of the data")
            if op + length > out_size:
                raise ReadError("PCD: LZF stream overruns the declared size")
            if ref + length <= op:
                out[op : op + length] = out[ref : ref + length]
            else:  # overlapping copy: the run repeats, byte by byte
                for k in range(length):
                    out[op + k] = out[ref + k]
            op += length
    if op != out_size:
        raise ReadError(f"PCD: LZF stream decoded to {op} bytes, expected {out_size}")
    return bytes(out)


def compress(data):
    """Encode ``data`` as an LZF stream."""
    src = bytes(data)
    n = len(src)
    out = bytearray()
    table = [-1] * (1 << _HASH_BITS)
    lit_start = 0
    ip = 0

    def flush(end):
        pos = lit_start
        while pos < end:
            run = min(_MAX_LITERAL, end - pos)
            out.append(run - 1)
            out.extend(src[pos : pos + run])
            pos += run

    while ip + 2 < n:
        h = ((src[ip] << 16) | (src[ip + 1] << 8) | src[ip + 2]) * 2654435761
        h = (h >> 8) & ((1 << _HASH_BITS) - 1)
        ref = table[h]
        table[h] = ip
        if (
            ref >= 0
            and ip - ref <= _MAX_OFFSET
            and src[ref] == src[ip]
            and src[ref + 1] == src[ip + 1]
            and src[ref + 2] == src[ip + 2]
        ):
            length = 3
            limit = min(_MAX_REF, n - ip)
            while length < limit and src[ref + length] == src[ip + length]:
                length += 1
            flush(ip)
            off = ip - ref - 1
            coded = length - 2
            if coded < 7:
                out.append((coded << 5) | (off >> 8))
            else:
                out.append((7 << 5) | (off >> 8))
                out.append(coded - 7)
            out.append(off & 0xFF)
            ip += length
            lit_start = ip
        else:
            ip += 1
    flush(n)
    return bytes(out)
