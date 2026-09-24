"""I/O for the LS-DYNA binary state database (the ``d3plot`` family), read-only.

The pure-Python twin of ``src/cpp/src/formats/lsdyna_d3plot.cpp``: both engines
read the same meshes.

A d3plot is word-addressed (4-byte words, or 8 in a double-precision run, in
either byte order): a control block of 64 words (plus ``EXTRA``), the geometry
(coordinates, then solid, thick-shell, beam and shell connectivity with a part
index as the last word), the user ids, optional part titles, then one state per
output time: the time, the global variables, the nodal variables, the element
variables of each family and the deletion flags. A run writes ``d3plot``,
``d3plot01``, ``d3plot02``... ; the family is one stream of states, which never
straddle a file. The word offsets and the state layout follow the LS-DYNA
database manual as lasso-python reads it. See ``doc/formats/lsdyna_d3plot.md``.
"""

import os
import re

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._mesh import CellBlock, Mesh
from .._regions import Region
from ..lsdyna._lsdyna import _collapse_solid

__all__ = ["read", "time_values", "family_files"]

TIME_KEY = "meshio:time"
_NAN = float("nan")
_FILETYPES = {1: "d3plot", 11: "d3eigv"}
_REFUSED_FILETYPES = {5: "d3part", 4: "intfor"}
_FEMZIP_NMMAT = 76_893_465
_EOF_MARKER = -999999.0


def _fail(message):
    raise ReadError(f"LS-DYNA d3plot: {message}")


def _digit(value, i):
    """The i-th decimal digit of ``value`` counted from the right (0 = units)."""
    return (abs(int(value)) // 10**i) % 10


# -- the family ------------------------------------------------------------------------


def family_files(filename):
    """The files of the family ``filename`` starts: itself, then ``<name>NN`` by number."""
    path = os.fspath(filename)
    folder = os.path.dirname(path) or "."
    base = os.path.basename(path)
    pattern = re.compile(re.escape(base) + r"([0-9]+)$")
    found = []
    try:
        entries = os.listdir(folder)
    except OSError:
        entries = []
    for entry in entries:
        m = pattern.match(entry)
        if m and os.path.isfile(os.path.join(folder, entry)):
            found.append((int(m.group(1)), entry))
    found.sort()
    return [path] + [os.path.join(os.path.dirname(path), e) for _, e in found]


def _continuation_base(path):
    """The base file when ``path`` is a numbered member of a family, else ``None``."""
    folder = os.path.dirname(path)
    m = re.match(r"^(.*?[^0-9])([0-9]+)$", os.path.basename(path))
    if not m:
        return None
    base = os.path.join(folder, m.group(1))
    return base if os.path.isfile(base) else None


# -- the control block -----------------------------------------------------------------


def _sniff(head):
    """(word size, byte order) of a d3plot control block, or ``None``."""
    for ws, order in ((4, "<"), (4, ">"), (8, "<"), (8, ">")):
        if len(head) < 64 * ws:
            continue
        itype = np.dtype(order + ("i4" if ws == 4 else "i8"))
        words = np.frombuffer(head[: 64 * ws], dtype=itype)
        filetype = int(words[11])
        if filetype > 1000:
            filetype -= 1000
        ndim = int(words[15])
        if filetype in (1, 4, 5, 11) and 2 <= ndim <= 9 and words[16] >= 0:
            return ws, order
    return None


_HEADER_WORDS = [
    (10, "runtime"),
    (11, "filetype"),
    (12, "source_version"),
    (15, "ndim"),
    (16, "numnp"),
    (17, "icode"),
    (18, "nglbv"),
    (19, "it"),
    (20, "iu"),
    (21, "iv"),
    (22, "ia"),
    (23, "nel8"),
    (24, "nummat8"),
    (25, "numds"),
    (26, "numst"),
    (27, "nv3d"),
    (28, "nel2"),
    (29, "nummat2"),
    (30, "nv1d"),
    (31, "nel4"),
    (32, "nummat4"),
    (33, "nv2d"),
    (34, "neiph"),
    (35, "neips"),
    (36, "maxint"),
    (37, "nmsph"),
    (38, "ngpsph"),
    (39, "narbs"),
    (40, "nelt"),
    (41, "nummatt"),
    (42, "nv3dt"),
    (43, "ioshl1"),
    (44, "ioshl2"),
    (45, "ioshl3"),
    (46, "ioshl4"),
    (47, "ialemat"),
    (48, "ncfdv1"),
    (49, "ncfdv2"),
    (50, "nadapt"),
    (51, "nmmat"),
    (52, "numfluid"),
    (53, "inn"),
    (54, "npefg"),
    (55, "nel48"),
    (56, "idtdt"),
    (57, "extra"),
]

_EXTRA_WORDS = [
    (64, "nel20"),
    (65, "nt3d"),
    (66, "nel27"),
    (67, "neipb"),
    (68, "nel21p"),
    (69, "nel15t"),
    (70, "soleng"),
    (71, "nel20t"),
    (72, "nel40p"),
    (73, "nel64"),
    (74, "quadr"),
    (75, "cubic"),
    (76, "tsheng"),
    (77, "nbranch"),
    (78, "penout"),
    (79, "engout"),
]


class _Words:
    """Word access into one file's bytes."""

    def __init__(self, data, ws, order):
        self.data = data
        self.ws = ws
        self.itype = np.dtype(order + ("i4" if ws == 4 else "i8"))
        self.ftype = np.dtype(order + ("f4" if ws == 4 else "f8"))

    def nwords(self):
        return len(self.data) // self.ws

    def ints(self, pos, n):
        if n <= 0:
            return np.zeros(0, dtype=np.int64)
        if (pos + n) * self.ws > len(self.data):
            _fail("the file is truncated")
        return np.frombuffer(self.data, self.itype, n, pos * self.ws).astype(np.int64)

    def floats(self, pos, n):
        if n <= 0:
            return np.zeros(0, dtype=np.float64)
        if (pos + n) * self.ws > len(self.data):
            _fail("the file is truncated")
        return np.frombuffer(self.data, self.ftype, n, pos * self.ws).astype(np.float64)

    def int(self, pos):
        return int(self.ints(pos, 1)[0])

    def float(self, pos):
        return float(self.floats(pos, 1)[0])

    def text(self, pos_bytes, nbytes):
        raw = self.data[pos_bytes : pos_bytes + nbytes]
        return raw.split(b"\0", 1)[0].decode("latin-1").strip(" \t\r\n\v\f")


class _Header:
    def __init__(self, w):
        h = {name: w.int(i) for i, name in _HEADER_WORDS}
        for i, name in _EXTRA_WORDS:
            h[name] = w.int(i) if h["extra"] > 0 and i < 64 + h["extra"] else 0
        self.raw = h
        self.title = w.text(0, 10 * w.ws)
        self.version = w.float(14)
        self.size = (64 + max(h["extra"], 0)) * w.ws  # bytes

        filetype = h["filetype"] - 1000 if h["filetype"] > 1000 else h["filetype"]
        if filetype in _REFUSED_FILETYPES:
            _fail(
                f"a {_REFUSED_FILETYPES[filetype]} file is not read (only d3plot/d3eigv)"
            )
        if filetype not in _FILETYPES:
            _fail(f"unknown file type {h['filetype']}")
        self.filetype = filetype
        if h["nmmat"] == _FEMZIP_NMMAT:
            _fail("the file is femzip-compressed; decompress it first")

        ndim = h["ndim"]
        self.has_material_type = ndim in (5, 7)
        self.has_rigid_road = ndim in (6, 7, 9)
        self.has_rigid_bodies = ndim in (8, 9)
        self.reduced_rigid_bodies = ndim == 9
        if ndim == 2:
            _fail("two-dimensional databases (NDIM = 2) are not read")
        self.nnodes = h["numnp"]
        self.nglbv = h["nglbv"]

        it = h["it"]
        self.node_mass_scaling = _digit(it, 1) == 1
        it0 = _digit(it, 0)
        self.node_temperature = it0 in (1, 2, 3)
        self.node_heat_flux = it0 in (2, 3)
        self.node_temperature_layers = it0 == 3
        self.node_displacement = h["iu"] != 0
        self.node_velocity = h["iv"] != 0
        self.node_acceleration = h["ia"] != 0

        self.nsolids = abs(h["nel8"])
        self.solid_extra_nodes = h["nel8"] < 0
        self.nbeams = h["nel2"]
        self.nshells = h["nel4"]
        self.ntshells = h["nelt"]
        self.nv3d = h["nv3d"]
        self.nv1d = h["nv1d"]
        self.nv2d = h["nv2d"]
        self.nv3dt = h["nv3dt"]
        self.neiph = h["neiph"]
        self.neips = h["neips"]
        self.neipb = h["neipb"]
        self.nt3d = h["nt3d"]

        maxint = h["maxint"]
        self.element_deletion = maxint <= -10000
        self.node_deletion = -10000 < maxint < 0
        self.layers = abs(maxint) - 10000 if maxint <= -10000 else abs(maxint)

        self.shell_stress = h["ioshl1"] == 1000
        self.solid_stress = h["ioshl1"] in (999, 1000)
        self.shell_pstrain = h["ioshl2"] == 1000
        self.solid_pstrain = h["ioshl2"] in (999, 1000)
        self.shell_forces = h["ioshl3"] == 1000
        self.shell_extra = h["ioshl4"] == 1000

        if h["ncfdv1"] != 0:
            _fail("CFD or multi-solver data (NCFDV1 != 0) is not read")
        if h["nadapt"] != 0:
            _fail("adaptive-remeshing databases (NADAPT != 0) are not read")
        self.nparts = h["nmmat"]
        self.nairbags = h["npefg"] % 1000 if 0 < h["npefg"] <= 10_000_000 else 0
        self.airbag_subver = h["npefg"] // 1000 if self.nairbags else 0
        self.nshells8 = h["nel48"]
        self.nsph = h["nmsph"]

        idtdt = h["idtdt"]
        self.node_temperature_gradient = _digit(idtdt, 0) == 1
        self.node_residual_forces = _digit(idtdt, 1) == 1
        self.plastic_strain_tensor = _digit(idtdt, 2) == 1
        self.thermal_strain_tensor = _digit(idtdt, 3) == 1
        if idtdt > 100:
            self.element_strain = _digit(idtdt, 4) == 1
        elif self.nv2d > 0:
            self.element_strain = (
                self.nv2d
                - self.layers
                * (6 * self.shell_stress + self.shell_pstrain + self.neips)
                - 8 * self.shell_forces
                - 4 * self.shell_extra
            ) > 1
        elif self.nv3dt > 0:
            self.element_strain = (
                self.nv3dt
                - self.layers
                * (6 * self.shell_stress + self.shell_pstrain + self.neips)
            ) > 1
        else:
            self.element_strain = False

        for key in ("nel20", "nel27", "nel21p", "nel15t", "nel20t", "nel40p", "nel64"):
            if h[key] > 0:
                _fail(f"higher-order solids ({key.upper()} = {h[key]}) are not read")

    @property
    def solid_layers(self):
        base = 6 * self.solid_stress + self.solid_pstrain + self.neiph
        return 8 if self.nv3d // max(base, 1) >= 8 else 1


# -- the geometry ----------------------------------------------------------------------


class _Geometry:
    """Everything between the control block and the first state."""

    def __init__(self, w, hdr):
        pos = hdr.size // w.ws
        ws = w.ws
        self.nrigid_shells = 0
        self.part_mattype = None
        if hdr.has_material_type:
            self.nrigid_shells = w.int(pos)
            nummat = w.int(pos + 1)
            if nummat != hdr.nparts:
                _fail(
                    f"the material type section lists {nummat} parts, the header {hdr.nparts}"
                )
            self.part_mattype = w.ints(pos + 2, hdr.nparts)
            pos += 2 + hdr.nparts
        if hdr.raw["ialemat"] > 0:
            pos += hdr.raw["ialemat"]
        self.nsph_vars = 0
        if hdr.nsph > 0:
            flags = w.ints(pos, 11)
            history = 0 if flags[0] == 10 else int(flags[10])
            self.nsph_vars = (
                int(flags[1:8].sum()) + abs(int(flags[8])) + int(flags[9]) + history + 1
            )
            pos += int(flags[0])
        self.airbag = None
        if hdr.nairbags:
            ngeom, nvar, npart, nstgeom = (int(v) for v in w.ints(pos, 4))
            pos += 4
            if hdr.airbag_subver == 4:
                pos += 1
            nvars = ngeom + nvar + nstgeom
            pos += nvars + 8 * nvars  # type codes, then 8-word names
            self.airbag = (ngeom, nvar, npart, nstgeom)

        # Coordinates and connectivity.
        n = hdr.nnodes
        self.coords = w.floats(pos, 3 * n).reshape(n, 3)
        pos += 3 * n
        self.solids = w.ints(pos, 9 * hdr.nsolids).reshape(hdr.nsolids, 9)
        pos += 9 * hdr.nsolids
        self.solid_extra = None
        if hdr.solid_extra_nodes:
            self.solid_extra = w.ints(pos, 2 * hdr.nsolids).reshape(hdr.nsolids, 2)
            pos += 2 * hdr.nsolids
        self.tshells = w.ints(pos, 9 * hdr.ntshells).reshape(hdr.ntshells, 9)
        pos += 9 * hdr.ntshells
        self.beams = w.ints(pos, 6 * hdr.nbeams).reshape(hdr.nbeams, 6)
        pos += 6 * hdr.nbeams
        self.shells = w.ints(pos, 5 * hdr.nshells).reshape(hdr.nshells, 5)
        pos += 5 * hdr.nshells

        # User ids.
        self.node_ids = np.arange(1, n + 1, dtype=np.int64)
        self.solid_ids = np.arange(1, hdr.nsolids + 1, dtype=np.int64)
        self.beam_ids = np.arange(1, hdr.nbeams + 1, dtype=np.int64)
        self.shell_ids = np.arange(1, hdr.nshells + 1, dtype=np.int64)
        self.tshell_ids = np.arange(1, hdr.ntshells + 1, dtype=np.int64)
        self.part_ids = None
        self.nrigid_bodies = 0
        narbs = hdr.raw["narbs"]
        if narbs > 0:
            start = pos
            head = w.ints(pos, 10)
            pos += 10
            nsort = int(head[0])
            counts = [
                int(v) for v in head[5:10]
            ]  # nodes, solids, beams, shells, tshells
            nparts_ids = hdr.nparts
            if nsort < 0:
                extra = w.ints(pos, 6)
                pos += 6
                nparts_ids = int(extra[3])
                self.nrigid_bodies = int(extra[4])
            ids = []
            for count in counts:
                ids.append(w.ints(pos, count))
                pos += count
            if len(ids[0]) == n:
                self.node_ids = ids[0]
            for attr, got, expected in (
                ("solid_ids", ids[1], hdr.nsolids),
                ("beam_ids", ids[2], hdr.nbeams),
                ("shell_ids", ids[3], hdr.nshells),
                ("tshell_ids", ids[4], hdr.ntshells),
            ):
                if len(got) == expected:
                    setattr(self, attr, got)
            if nsort < 0 and nparts_ids == hdr.nparts:
                self.part_ids = w.ints(pos, hdr.nparts)
            pos = start + narbs

        # Rigid bodies (skipped: their motion is not read).
        if hdr.has_rigid_bodies:
            nrigid = w.int(pos)
            pos += 1
            for _ in range(nrigid):
                numnodr = w.int(pos + 1)
                pos += 2 + numnodr
                numnoda = w.int(pos)
                pos += 1 + numnoda
            self.nrigid_bodies_motion = nrigid
        else:
            self.nrigid_bodies_motion = 0
        if hdr.nsph > 0:
            pos += 2 * hdr.nsph
        if self.airbag is not None:
            ngeom = self.airbag[0]
            pos += hdr.nairbags * ngeom
        self.nroads = 0
        if hdr.has_rigid_road:
            nnode, nseg, nsurf = (int(v) for v in w.ints(pos, 3))
            pos += 4 + nnode + 3 * nnode
            for _ in range(nsurf):
                surf_nseg = w.int(pos + 1)
                pos += 2 + 4 * surf_nseg
            self.nroads = nsurf

        # Extra connectivity. The ten-node solids' two extra nodes follow the
        # 8-node connectivity (database manual); lasso-python also writes and
        # reads a second, identical copy here, which is skipped when present.
        self.tet_extra = self.solid_extra
        if hdr.solid_extra_nodes and hdr.nsolids > 0:
            n2 = 2 * hdr.nsolids
            if (pos + n2) * ws <= len(w.data) and np.array_equal(
                w.ints(pos, n2), self.solid_extra.reshape(-1)
            ):
                pos += n2
        self.shell8 = None
        if hdr.nshells8 > 0:
            self.shell8 = w.ints(pos, 5 * hdr.nshells8).reshape(hdr.nshells8, 5)
            pos += 5 * hdr.nshells8

        # Part titles (and the end-of-geometry marker in front of them).
        self.part_titles = {}
        if pos < w.nwords() and w.float(pos) == _EOF_MARKER:
            pos += 1
            while pos < w.nwords():
                ntype = w.int(pos)
                if ntype == 90000:
                    pos += 1 + 72 // ws
                elif ntype in (90001, 90002, 90020):
                    count = w.int(pos + 1)
                    pos += 2
                    entry = ws + 72
                    for k in range(count):
                        at = pos * ws + k * entry
                        if ntype == 90001:
                            pid = w.int((at) // ws)
                            self.part_titles.setdefault(pid, w.text(at + ws, 72))
                    pos += (count * entry) // ws
                else:
                    break
                if pos < w.nwords() and w.float(pos) == _EOF_MARKER:
                    pos += 1
        self.end = pos  # words


# -- the state layout ------------------------------------------------------------------


def _node_vars(hdr):
    """The nodal variables of a state in file order: (name, components)."""
    out = []
    if hdr.node_displacement:
        out.append(("coordinates", 3))
    if hdr.node_temperature:
        out.append(("temperature", 3 if hdr.node_temperature_layers else 1))
    if hdr.node_heat_flux:
        out.append(("heat_flux", 3))
    if hdr.node_mass_scaling:
        out.append(("mass_scaling", 1))
    if hdr.node_temperature_gradient:
        out.append(("temperature_gradient", 1))
    if hdr.node_residual_forces:
        out.append(("residual_forces", 3))
        out.append(("residual_moments", 3))
    if hdr.node_velocity:
        out.append(("velocity", 3))
    if hdr.node_acceleration:
        out.append(("acceleration", 3))
    return out


def _state_words(hdr, geo):
    n = 1 + hdr.nglbv
    n += sum(c for _, c in _node_vars(hdr)) * hdr.nnodes
    n += hdr.nt3d * hdr.nsolids
    n += hdr.nsolids * hdr.nv3d
    n += hdr.ntshells * hdr.nv3dt
    n += hdr.nbeams * hdr.nv1d
    n += (hdr.nshells - geo.nrigid_shells) * hdr.nv2d
    n += hdr.nsph * geo.nsph_vars
    if hdr.node_deletion:
        n += hdr.nnodes
    elif hdr.element_deletion:
        n += hdr.nbeams + hdr.nshells + hdr.nsolids + hdr.ntshells
    if geo.airbag is not None:
        _, nvar, npart, nstgeom = geo.airbag
        n += hdr.nairbags * nstgeom + npart * nvar
    n += geo.nroads * 6
    if hdr.has_rigid_bodies:
        n += geo.nrigid_bodies_motion * (12 if hdr.reduced_rigid_bodies else 24)
    return n


def _last_nonzero_byte(path, size):
    """One past the last non-zero byte of ``path`` (0 for an all-zero file)."""
    chunk = 1 << 16
    with open(path, "rb") as f:
        end = size
        while end > 0:
            start = max(0, end - chunk)
            f.seek(start)
            block = np.frombuffer(f.read(end - start), dtype=np.uint8)
            nz = np.flatnonzero(block)
            if nz.size:
                return start + int(nz[-1]) + 1
            end = start
    return 0


class _File:
    """A d3plot family: control block, geometry and where each state lies."""

    def __init__(self, filename):
        path = os.fspath(filename)
        if not os.path.isfile(path):
            _fail(f"'{path}' does not exist")
        base = _continuation_base(path)
        if base is not None and _sniff_file(base):
            _fail(
                f"'{os.path.basename(path)}' continues the family of '{os.path.basename(base)}'"
                "; open the base file"
            )
        with open(path, "rb") as f:
            data = f.read()
        if path.lower().endswith(".fz"):
            _fail("the file is femzip-compressed; decompress it first")
        sniffed = _sniff(data)
        if sniffed is None:
            _fail(f"'{path}' is not a d3plot file (no plausible control block)")
        self.ws, self.order = sniffed
        w = _Words(data, *sniffed)
        self.words = w
        self.hdr = hdr = _Header(w)
        self.geo = geo = _Geometry(w, hdr)
        self.state_words = _state_words(hdr, geo)
        if self.state_words <= 0:
            _fail("the header describes an empty state")
        state_bytes = self.state_words * self.ws

        # (path, byte offset) of every state, family order.
        self.states = []
        for k, member in enumerate(family_files(path)):
            size = os.path.getsize(member)
            start = geo.end * self.ws if k == 0 else 0
            # the last non-zero *word*: a big-endian value can end in zero bytes
            last = -(-_last_nonzero_byte(member, size) // self.ws) * self.ws
            count = max(0, (last - start) // state_bytes)
            if hdr.filetype == 11 and k > 0 and count == 0 and state_bytes <= size:
                count = 1
            self.states.extend((member, start + i * state_bytes) for i in range(count))

    def state(self, index):
        member, offset = self.states[index]
        with open(member, "rb") as f:
            f.seek(offset)
            data = f.read(self.state_words * self.ws)
        if len(data) < self.state_words * self.ws:
            _fail(f"state {index} is truncated")
        return np.frombuffer(data, self.words.ftype).astype(np.float64)

    def times(self):
        out = []
        for member, offset in self.states:
            with open(member, "rb") as f:
                f.seek(offset)
                raw = f.read(self.ws)
            out.append(float(np.frombuffer(raw, self.words.ftype)[0]))
        return out


def _sniff_file(path):
    try:
        with open(path, "rb") as f:
            return _sniff(f.read(64 * 8)) is not None
    except OSError:
        return False


def is_d3plot(head):
    """Whether ``head`` (the first bytes of a file) starts a d3plot control block.

    Stricter than the reader's own sniff, since it runs on files of any format:
    a plausible file type and dimension, non-negative counts, a version of 0
    (lasso-python writes none) or between 900 and 100000, and a printable title.
    """
    sniffed = _sniff(head)
    if sniffed is None:
        return False
    ws, order = sniffed
    w = _Words(head, ws, order)
    version = w.float(14)
    if not (version == 0.0 or 900.0 <= version < 100000.0):
        return False
    if any(w.int(i) < 0 for i in (16, 18, 28, 31, 40)):
        return False
    return all(b == 0 or 32 <= b < 127 for b in head[:40])


# -- the mesh --------------------------------------------------------------------------

# families in block order: (name, dimension)
_SOLID, _TSHELL, _BEAM, _SHELL = range(4)


def _build_cells(hdr, geo):
    """Cell blocks, one per (family, cell type) in first-appearance order.

    Returns the blocks, per family the (block, row) of each element, and per
    block the family, element indices and part indices.
    """
    blocks = []  # [type, rows, family, elems, parts]
    where = [None] * 4

    def add(family, elements, part_word, make):
        index = {}
        pos = np.empty((len(elements), 2), dtype=np.int64)
        for e, row in enumerate(elements):
            cell_type, nodes = make(e, row)
            key = cell_type
            if key not in index:
                index[key] = len(blocks)
                blocks.append([cell_type, [], family, [], []])
            b = blocks[index[key]]
            pos[e] = (index[key], len(b[1]))
            b[1].append(nodes)
            b[3].append(e)
            b[4].append(int(row[part_word]) - 1)
        where[family] = pos

    def solid(e, row):
        nodes = [int(v) - 1 for v in row[:8]]
        if (
            geo.tet_extra is not None
            and geo.tet_extra[e, 0] > 0
            and geo.tet_extra[e, 1] > 0
        ):
            return "tetra10", nodes + [int(v) - 1 for v in geo.tet_extra[e]]
        cell_type, kept = _collapse_solid(nodes)
        return cell_type, kept

    def beam(e, row):
        return "line", [int(row[0]) - 1, int(row[1]) - 1]

    shell8 = {}
    if geo.shell8 is not None:
        for row in geo.shell8:
            shell8[int(row[0]) - 1] = [int(v) - 1 for v in row[1:5]]

    def shell(e, row):
        nodes = [int(v) - 1 for v in row[:4]]
        if e in shell8:
            return "quad8", nodes + shell8[e]
        if nodes[3] == nodes[2] or nodes[3] < 0:
            return "triangle", nodes[:3]
        return "quad", nodes

    add(_SOLID, geo.solids, 8, solid)
    add(_TSHELL, geo.tshells, 8, solid)
    add(_BEAM, geo.beams, 5, beam)
    add(_SHELL, geo.shells, 4, shell)
    return blocks, where


def _part_user_ids(hdr, geo):
    if geo.part_ids is not None:
        return [int(v) for v in geo.part_ids]
    titles = list(geo.part_titles)
    if len(titles) == hdr.nparts:
        return titles
    return list(range(1, hdr.nparts + 1))


_FAMILY_DIM = {_SOLID: 3, _TSHELL: 3, _BEAM: 1, _SHELL: 2}


def _dim(cell_type):
    return {"line": 1, "triangle": 2, "quad": 2, "quad8": 2}.get(cell_type, 3)


def _build_mesh(f):
    hdr, geo = f.hdr, f.geo
    nnodes = hdr.nnodes
    points = geo.coords.copy()
    blocks, where = _build_cells(hdr, geo)
    for b in blocks:
        rows = np.array(b[1], dtype=np.int64).reshape(len(b[1]), -1)
        if rows.size and (rows.min() < 0 or rows.max() >= nnodes):
            _fail("an element references a node outside the node table")
    cells = [CellBlock(b[0], np.array(b[1], dtype=np.int64)) for b in blocks]
    mesh = Mesh(points, cells)
    mesh.point_data["lsdyna:nid"] = geo.node_ids.astype(np.int64)
    family_ids = {
        _SOLID: geo.solid_ids,
        _TSHELL: geo.tshell_ids,
        _BEAM: geo.beam_ids,
        _SHELL: geo.shell_ids,
    }
    part_ids = _part_user_ids(hdr, geo)
    if blocks:
        mesh.cell_data["lsdyna:eid"] = [
            family_ids[b[2]][np.array(b[3], dtype=np.int64)].astype(np.int64)
            for b in blocks
        ]

        def user_part(p):
            return part_ids[p] if 0 <= p < len(part_ids) else p + 1

        mesh.cell_data["lsdyna:part"] = [
            np.array([user_part(p) for p in b[4]], dtype=np.int64) for b in blocks
        ]
        by_part = {}
        base = 0
        for b in blocks:
            d = _dim(b[0])
            for i, p in enumerate(b[4]):
                entry = by_part.setdefault(user_part(p), [-1, []])
                entry[0] = max(entry[0], d)
                entry[1].append(base + i)
            base += len(b[1])
        regions = []
        for pid, (d, members) in by_part.items():
            title = geo.part_titles.get(pid, "")
            name = title if title else f"Part {pid}"
            regions.append(
                Region(name, "cell", np.array(members, dtype=np.int64), dim=d, tag=pid)
            )
        regions.sort(key=lambda r: (1, r.name, r.dim, r.tag))
        mesh.regions = regions
    return mesh, blocks, where


# -- one state -------------------------------------------------------------------------


class _CellArrays:
    """Named per-cell arrays, widths padded with NaN across blocks."""

    def __init__(self, blocks):
        self.blocks = blocks
        self.values = {}  # name -> (points, width, per-family (elements, rows) list)
        self.order = []

    def add(self, name, family, rows, points, width):
        """``rows``: (nelements, points * width) values of ``family`` in element order."""
        if name not in self.values:
            self.values[name] = [points, width, []]
            self.order.append(name)
        entry = self.values[name]
        entry[0] = max(entry[0], points)
        entry[1] = max(entry[1], width)
        entry[2].append((family, rows, points, width))

    def emit(self, mesh, where, fmt="lsdyna_d3plot"):
        for name in self.order:
            pts, width, parts = self.values[name]
            cols = pts * width
            out = []
            for b in self.blocks:
                a = np.full((len(b[1]), cols), _NAN)
                out.append(a)
            for family, rows, p, wdt in parts:
                pos = where[family]
                if pos is None or len(pos) == 0:
                    continue
                rows = rows.reshape(len(rows), p, wdt)
                padded = np.full((len(rows), pts, width), _NAN)
                padded[:, :p, :wdt] = rows
                padded = padded.reshape(len(rows), cols)
                for e in range(len(rows)):
                    bi, ri = pos[e]
                    out[bi][ri] = padded[e]
            if cols == 1:
                out = [a.reshape(-1) for a in out]
            mesh.cell_data[name] = out
            if pts > 1:
                mesh.field_data[f"{fmt}:layout:{name}"] = np.array(
                    [pts, width], dtype=np.int64
                )


def _read_state(f, mesh, blocks, where, index, wanted):
    hdr, geo = f.hdr, f.geo
    s = f.state(index)
    nn = hdr.nnodes
    k = 0
    time = s[k]
    k += 1
    globals_ = s[k : k + hdr.nglbv]
    k += hdr.nglbv

    def want(name):
        return wanted is None or name in wanted

    # globals: kinetic, internal, total energy, velocity, then per part (internal,
    # kinetic energy, velocity, mass, hourglass energy), then rigid walls.
    g = 0
    for name, width in (
        ("global_kinetic_energy", 1),
        ("global_internal_energy", 1),
        ("global_total_energy", 1),
        ("global_velocity", 3),
    ):
        if g + width <= hdr.nglbv:
            if want(name):
                mesh.field_data[name] = globals_[g : g + width].copy()
            g += width
    nparts = hdr.raw["nummat8"] + hdr.raw["nummat2"] + hdr.raw["nummat4"]
    nparts += hdr.raw["nummatt"] + geo.nrigid_bodies
    for name, width in (
        ("part_internal_energy", 1),
        ("part_kinetic_energy", 1),
        ("part_velocity", 3),
        ("part_mass", 1),
        ("part_hourglass_energy", 1),
    ):
        if g + width * nparts <= hdr.nglbv:
            if want(name):
                v = globals_[g : g + width * nparts]
                mesh.field_data[name] = v.reshape(nparts, 3) if width == 3 else v.copy()
            g += width * nparts

    # nodes
    for name, comps in _node_vars(hdr):
        v = (
            s[k : k + comps * nn].reshape(nn, comps)
            if comps > 1
            else s[k : k + nn].copy()
        )
        k += comps * nn
        if name == "coordinates":
            if want("displacement"):
                mesh.point_data["displacement"] = v - geo.coords
        elif want(name):
            mesh.point_data[name] = v

    cells = _CellArrays(blocks)

    def put(name, family, rows, points, width):
        if want(name):
            cells.add(name, family, rows, points, width)

    # solid thermal
    if hdr.nt3d > 0:
        n = hdr.nsolids
        put(
            "thermal_variables",
            _SOLID,
            s[k : k + n * hdr.nt3d].reshape(n, -1),
            1,
            hdr.nt3d,
        )
        k += n * hdr.nt3d

    # solids
    if hdr.nsolids > 0 and hdr.nv3d > 0:
        n, nv = hdr.nsolids, hdr.nv3d
        layers = hdr.solid_layers
        data = s[k : k + n * nv].reshape(n, layers, nv // layers)
        i = 0
        if hdr.solid_stress:
            put("stress", _SOLID, data[:, :, i : i + 6].reshape(n, -1), layers, 6)
            i += 6
        if hdr.solid_pstrain:
            put("effective_plastic_strain", _SOLID, data[:, :, i], layers, 1)
            i += 1
        history = data[:, :, i : i + hdr.neiph]
        nstrain = 6 * hdr.element_strain
        if nstrain and history.shape[2] >= nstrain:
            put("strain", _SOLID, history[:, :, -nstrain:].reshape(n, -1), layers, 6)
            history = history[:, :, :-nstrain]
        if hdr.plastic_strain_tensor and history.shape[2] >= 6:
            put(
                "plastic_strain_tensor",
                _SOLID,
                history[:, :, :6].reshape(n, -1),
                layers,
                6,
            )
            history = history[:, :, 6:]
        if hdr.thermal_strain_tensor and history.shape[2] >= 6:
            put(
                "thermal_strain_tensor",
                _SOLID,
                history[:, :, :6].reshape(n, -1),
                layers,
                6,
            )
            history = history[:, :, 6:]
        if history.shape[2]:
            put(
                "history_variables",
                _SOLID,
                history.reshape(n, -1),
                layers,
                history.shape[2],
            )
        k += n * nv

    # thick shells
    if hdr.ntshells > 0 and hdr.nv3dt > 0:
        n, nv, nl = hdr.ntshells, hdr.nv3dt, hdr.layers
        nh = hdr.neips
        nlayer = nl * (6 * hdr.shell_stress + hdr.shell_pstrain + nh)
        data = s[k : k + n * nv].reshape(n, nv)
        layer = (
            data[:, :nlayer].reshape(n, nl, -1) if nl else data[:, :0].reshape(n, 0, 0)
        )
        i = 0
        if hdr.shell_stress:
            put("stress", _TSHELL, layer[:, :, i : i + 6].reshape(n, -1), nl, 6)
            i += 6
        if hdr.shell_pstrain:
            put("effective_plastic_strain", _TSHELL, layer[:, :, i], nl, 1)
            i += 1
        if nh:
            put(
                "history_variables",
                _TSHELL,
                layer[:, :, i : i + nh].reshape(n, -1),
                nl,
                nh,
            )
        if hdr.element_strain:
            strain = data[:, nlayer : nlayer + 12]
            put("strain_inner", _TSHELL, strain[:, :6], 1, 6)
            put("strain_outer", _TSHELL, strain[:, 6:12], 1, 6)
        k += n * nv

    # beams
    if hdr.nbeams > 0 and hdr.nv1d > 0:
        n, nv, nh = hdr.nbeams, hdr.nv1d, hdr.neipb
        nl = int((-3 * nh + nv - 6) / (nh + 5))
        data = s[k : k + n * nv].reshape(n, nv)
        put("beam_axial_force", _BEAM, data[:, 0], 1, 1)
        put("beam_shear_force", _BEAM, data[:, 1:3], 1, 2)
        put("beam_bending_moment", _BEAM, data[:, 3:5], 1, 2)
        put("beam_torsion_moment", _BEAM, data[:, 5], 1, 1)
        if nl > 0:
            layer = data[:, 6 : 6 + 5 * nl].reshape(n, nl, 5)
            put("beam_axial_stress", _BEAM, layer[:, :, 0], nl, 1)
            put("beam_shear_stress", _BEAM, layer[:, :, 1:3].reshape(n, -1), nl, 2)
            put("effective_plastic_strain", _BEAM, layer[:, :, 3], nl, 1)
            put("beam_axial_strain", _BEAM, layer[:, :, 4], nl, 1)
        if nh:
            hist = data[:, 6 + 5 * nl :]
            put("history_variables", _BEAM, hist, 3 + nl, nh)
        k += n * nv

    # shells (rigid ones have no state data)
    nreduced = hdr.nshells - geo.nrigid_shells
    if nreduced > 0 and hdr.nv2d > 0:
        nv, nl, nh = hdr.nv2d, hdr.layers, hdr.neips
        nlayer = nl * (6 * hdr.shell_stress + hdr.shell_pstrain + nh)
        data = s[k : k + nreduced * nv].reshape(nreduced, nv)
        k += nreduced * nv
        if geo.nrigid_shells:
            deformable = geo.part_mattype[geo.shells[:, 4] - 1] != 20
            full = np.full((hdr.nshells, nv), _NAN)
            full[deformable] = data
            data = full
        n = hdr.nshells
        layer = (
            data[:, :nlayer].reshape(n, nl, -1) if nl else data[:, :0].reshape(n, 0, 0)
        )
        i = 0
        if hdr.shell_stress:
            put("stress", _SHELL, layer[:, :, i : i + 6].reshape(n, -1), nl, 6)
            i += 6
        if hdr.shell_pstrain:
            put("effective_plastic_strain", _SHELL, layer[:, :, i], nl, 1)
            i += 1
        if nh:
            put(
                "history_variables",
                _SHELL,
                layer[:, :, i : i + nh].reshape(n, -1),
                nl,
                nh,
            )
        rest = data[:, nlayer:]
        j = 0
        if hdr.shell_forces:
            put("shell_bending_moment", _SHELL, rest[:, 0:3], 1, 3)
            put("shell_shear_force", _SHELL, rest[:, 3:5], 1, 2)
            put("shell_normal_force", _SHELL, rest[:, 5:8], 1, 3)
            j = 8
        if hdr.shell_extra:
            put("thickness", _SHELL, rest[:, j], 1, 1)
            put("shell_element_variables", _SHELL, rest[:, j + 1 : j + 3], 1, 2)
            j += 3
        if hdr.element_strain:
            put("strain_inner", _SHELL, rest[:, j : j + 6], 1, 6)
            put("strain_outer", _SHELL, rest[:, j + 6 : j + 12], 1, 6)
            j += 12
        if hdr.shell_extra:
            put("internal_energy", _SHELL, rest[:, j], 1, 1)
            j += 1
        if hdr.plastic_strain_tensor:
            put(
                "plastic_strain_tensor",
                _SHELL,
                rest[:, j : j + 6 * nl],
                nl,
                6,
            )
            j += 6 * nl
        if hdr.thermal_strain_tensor:
            put("thermal_strain_tensor", _SHELL, rest[:, j : j + 6], 1, 6)
            j += 6

    # SPH particles are skipped
    k += hdr.nsph * geo.nsph_vars

    # deletion
    if hdr.node_deletion:
        if want("lsdyna:alive"):
            mesh.point_data["lsdyna:alive"] = (s[k : k + nn] != 0).astype(np.int8)
        k += nn
    elif hdr.element_deletion:
        alive = {}
        for family, n in (
            (_SOLID, hdr.nsolids),
            (_TSHELL, hdr.ntshells),
            (_SHELL, hdr.nshells),
            (_BEAM, hdr.nbeams),
        ):
            alive[family] = s[k : k + n] != 0
            k += n
        if want("lsdyna:alive") and blocks:
            out = []
            for b in blocks:
                out.append(alive[b[2]][np.array(b[3], dtype=np.int64)].astype(np.int8))
            mesh.cell_data["lsdyna:alive"] = out

    cells.emit(mesh, where)
    return time


# -- entry points ----------------------------------------------------------------------


def time_values(filename):
    """The time of every state of the family."""
    return _File(filename).times()


def read(filename, points_only=False, arrays=None, time_step=0):
    f = _File(filename)
    mesh, blocks, where = _build_mesh(f)
    n = len(f.states)
    if n == 0:
        if time_step not in (0, -1):
            raise ReadError(
                f"time step {time_step} is out of range: the file has no states"
            )
        return mesh
    index = time_step + n if time_step < 0 else time_step
    if not 0 <= index < n:
        raise ReadError(
            f"time step {time_step} is out of range: the file has {n} step(s)"
        )
    mesh.time_values = f.times()
    mesh.field_data[TIME_KEY] = np.array([mesh.time_values[index]], dtype=np.float64)
    mesh.field_data["lsdyna:state"] = np.array([index], dtype=np.int64)
    if points_only:
        return mesh
    wanted = None if arrays is None else set(arrays)
    _read_state(f, mesh, blocks, where, index, wanted)
    if f.hdr.airbag_subver or f.hdr.nsph:
        warn("LS-DYNA d3plot: airbag particle and SPH data are skipped")
    return mesh
