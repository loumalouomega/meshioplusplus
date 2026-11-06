"""Ragged cell blocks that carry their CSR storage (roadmap 5.6).

The core hands a ragged block (jagged polygons, polyhedra) to Python as a list
of 1-D int64 views, one per polygon or face, of one flat buffer. Giving it back
used to visit every row from Python. :class:`CsrRows` is that list, unchanged
for every reader of it, which also keeps the flat buffer and its offsets, so
the core takes the block back with one copy of each and the Python layers can
vectorise over the arrays.

The arrays are dropped as soon as the list changes shape: every list mutator
(``rows[i] = ...``, ``append``, ``del``, ``sort``, ...) on the block, or on one
of a polyhedron's face lists (:class:`CsrFaces`), forgets them, and the block
is then read row by row as any other list. Writing into a row (``rows[i][0] =
7``) writes into the flat buffer the row views, so the two stay one.
"""

import numpy as np

__all__ = ["CsrFaces", "CsrRows", "csr_of"]


def _dropping(name):
    base = getattr(list, name)

    def method(self, *args, **kwargs):
        self._state[0] = None
        return base(self, *args, **kwargs)

    method.__name__ = name
    method.__doc__ = base.__doc__
    return method


class _Tracked(list):
    # ``_state`` is a one-element list shared by a block and its face lists:
    # ``[(flat, offsets, face_offsets)]``, or ``[None]`` once any of them has
    # changed shape. A shared cell rather than a reference to the parent, so
    # the face lists hold no cycle and the block is freed by reference count.
    __slots__ = ("_state",)

    def __init__(self, iterable=()):
        super().__init__(iterable)
        self._state = [None]

    __setitem__ = _dropping("__setitem__")
    __delitem__ = _dropping("__delitem__")
    __iadd__ = _dropping("__iadd__")
    __imul__ = _dropping("__imul__")
    append = _dropping("append")
    extend = _dropping("extend")
    insert = _dropping("insert")
    pop = _dropping("pop")
    remove = _dropping("remove")
    clear = _dropping("clear")
    sort = _dropping("sort")
    reverse = _dropping("reverse")


class CsrFaces(_Tracked):
    """One polyhedron's faces in a :class:`CsrRows` block; changing this list
    drops the block's CSR arrays."""

    __slots__ = ()

    def __reduce__(self):
        return (CsrFaces, (list(self),))


class CsrRows(_Tracked):
    """A ragged block's rows, a ``list``, with the CSR arrays behind them.

    ``csr()`` returns ``(flat, offsets)`` for polygons and ``(flat, offsets,
    face_offsets)`` for polyhedra, all int64: row ``r`` is
    ``flat[offsets[r]:offsets[r + 1]]``, and polyhedron ``c`` has the faces
    (rows) ``face_offsets[c]`` to ``face_offsets[c + 1]``. It returns ``None``
    once the list has changed shape.
    """

    __slots__ = ()

    @classmethod
    def from_csr(cls, flat, offsets, face_offsets=None):
        """The block over ``flat``: each row is a view of it."""
        flat = np.ascontiguousarray(flat, dtype=np.int64)
        offsets = np.ascontiguousarray(offsets, dtype=np.int64)
        rows = np.split(flat, offsets[1:-1]) if len(offsets) > 1 else []
        state = [(flat, offsets, None)]
        if face_offsets is None:
            out = cls(rows)
        else:
            face_offsets = np.ascontiguousarray(face_offsets, dtype=np.int64)
            state = [(flat, offsets, face_offsets)]
            out = cls()
            list.extend(
                out,
                (
                    CsrFaces(rows[face_offsets[c] : face_offsets[c + 1]])
                    for c in range(len(face_offsets) - 1)
                ),
            )
            for cell in out:
                cell._state = state
        out._state = state
        return out

    def csr(self):
        """The CSR arrays behind the rows, or ``None`` (see the class)."""
        return self._state[0]

    def __reduce__(self):
        # Rebuilt from the arrays, so a copy's rows view the copy's buffer
        # (``copy.deepcopy`` copies the arrays first; ``copy.copy`` shares).
        csr = self._state[0]
        if csr is None:
            return (CsrRows, (list(self),))
        return (CsrRows.from_csr, csr)


def csr_of(data):
    """The CSR arrays a ragged block's ``data`` still carries, or ``None``."""
    return data.csr() if isinstance(data, CsrRows) else None
