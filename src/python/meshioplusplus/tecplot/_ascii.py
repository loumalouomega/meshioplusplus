"""
ASCII Tecplot scanner: VARIABLES and every ZONE header, each zone's data
located by counting its own tokens (the Python twin of ``tecplot_scan_zones``
and ``tecplot_decode_zone`` in ``tecplot.cpp``).
"""

import numpy as np

from .._exceptions import ReadError
from .._files import open_file
from ._zones import Zone, face_map

_ET = {
    "LINESEG": "FELINESEG",
    "TRIANGLE": "FETRIANGLE",
    "QUADRILATERAL": "FEQUADRILATERAL",
    "TETRAHEDRON": "FETETRAHEDRON",
    "BRICK": "FEBRICK",
}


def _is_float(tok):
    try:
        float(tok)
    except ValueError:
        return False
    return True


def header_tokens(s):
    """``tecplot_header_tokens``: quoted and parenthesised spans stay one token,
    ``=`` is its own token, everything else splits on blanks and commas."""
    out = []
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c in " \t,":
            i += 1
        elif c == "=":
            out.append("=")
            i += 1
        elif c in '"(':
            j = s.find('"' if c == '"' else ")", i + 1)
            if j < 0:
                j = n - 1
            out.append(s[i : j + 1])
            i = j + 1
        else:
            j = i
            while j < n and s[j] not in ' \t,="(':
                j += 1
            out.append(s[i:j])
            i = j
    return out


def _unquote(tok):
    if len(tok) >= 2 and tok[0] == '"' and tok[-1] == '"':
        return tok[1:-1]
    return tok


def split_entries(inner):
    """Splits on commas outside ``[...]``."""
    out, cur, depth = [], "", 0
    for c in inner:
        if c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
        if c == "," and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += c
    if cur:
        out.append(cur)
    return out


def parse_ranges(spec):
    """``[a,b-c]`` (1-based) -> 0-based indices."""
    spec = spec.strip()
    if spec.startswith("["):
        spec = spec[1:]
    if spec.endswith("]"):
        spec = spec[:-1]
    out = []
    for part in split_entries(spec):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a) - 1, int(b)))
        else:
            out.append(int(part) - 1)
    return out


def _variables(line):
    rhs = line[line.find("=") + 1 :]
    out = []
    p = 0
    while p < len(rhs):
        c = rhs[p]
        if c == '"':
            q = rhs.find('"', p + 1)
            if q < 0:
                q = len(rhs)
            out.append(rhs[p + 1 : q])
            p = q + 1
        elif c.isspace() or c == ",":
            p += 1
        else:
            q = p
            while q < len(rhs) and not rhs[q].isspace() and rhs[q] != ",":
                q += 1
            out.append(rhs[p:q])
            p = q
    return out


def _zone_kind(z, fields, nvar):
    """Zone type, packing and cell-centred variables from the header fields."""
    f = fields.get("F", "").upper()
    if f:
        if f in ("FEPOINT", "FEBLOCK"):
            et = fields.get("ET", "").upper()
            z.type_name = _ET.get(et, et)
            z.block = f == "FEBLOCK"
        else:
            z.type_name = "ORDERED"
            z.block = f == "BLOCK"
    else:
        z.type_name = fields.get("ZONETYPE", "ORDERED").upper()
        z.block = fields.get("DATAPACKING", "BLOCK").upper() == "BLOCK"
    z.ordered = z.type_name == "ORDERED"
    z.cell_centered = [0] * nvar
    if not z.block:
        return
    if "NV" in fields and not z.ordered:
        for k in range(int(fields["NV"]), nvar):
            z.cell_centered[k] = 1
    elif z.varloc:
        for entry in split_entries(z.varloc[1:-1]):
            if "=" not in entry:
                continue
            rng, loc = entry.split("=", 1)
            if loc.strip().upper() != "CELLCENTERED":
                continue
            for idx in parse_ranges(rng):
                if 0 <= idx < nvar:
                    z.cell_centered[idx] = 1


def _data_tokens(z, nvar):
    if not z.block:
        return z.num_nodes * nvar
    return sum(z.data_length(v) for v in range(nvar) if z.owns(v))


def _face_map_tokens(z):
    """How many integers a face-based zone's face map holds: node counts
    (polyhedra), face nodes, left and right elements, boundary connections."""
    n = (z.num_faces if z.is_polyhedron else 0) + z.total_face_nodes + 2 * z.num_faces
    if z.num_boundary_faces > 0:
        n += z.num_boundary_faces + 2 * z.num_boundary_conns
    return n


def _tokens(line):
    return line.split()


class AsciiSource:
    def __init__(self, lines, zones, nvar):
        self.lines = lines
        self.zones = zones
        self.nvar = nvar

    def own_data(self, idx):
        z = self.zones[idx]
        nvar = self.nvar
        want = _data_tokens(z, nvar)
        flat = []
        li = z.data_start
        while len(flat) < want and li < len(self.lines):
            flat.extend(float(t) for t in _tokens(self.lines[li]))
            li += 1
        if len(flat) < want:
            raise ReadError(
                f"Tecplot: zone {idx + 1} has fewer values than its header announces"
            )
        flat = np.array(flat[:want]) if want else np.zeros(0)
        cols = {}
        if z.block:
            off = 0
            for v in range(nvar):
                if not z.owns(v):
                    continue
                n = z.data_length(v)
                cols[v] = flat[off : off + n]
                off += n
        else:
            table = flat.reshape(z.num_nodes, nvar)
            for v in range(nvar):
                if z.owns(v):
                    cols[v] = table[:, v].copy()
        conn = None
        if not z.ordered and z.conn_share < 0 and z.is_poly:
            want = _face_map_tokens(z)
            ints = []
            while len(ints) < want and li < len(self.lines):
                ints.extend(int(t) for t in _tokens(self.lines[li]))
                li += 1
            if len(ints) < want:
                raise ReadError(f"Tecplot: zone {idx + 1} face map is truncated")
            pos = 0

            def take(n):
                nonlocal pos
                pos += n
                return ints[pos - n : pos]

            counts = take(z.num_faces) if z.is_polyhedron else None
            nodes = take(z.total_face_nodes)
            left = take(z.num_faces)
            right = take(z.num_faces)
            return cols, face_map(z, idx, counts, nodes, left, right, True)
        if not z.ordered and z.conn_share < 0:
            npc = {
                "FELINESEG": 2,
                "FETRIANGLE": 3,
                "FEQUADRILATERAL": 4,
                "FETETRAHEDRON": 4,
            }.get(z.type_name, 8)
            rows = []
            for _ in range(z.num_cells):
                if li >= len(self.lines):
                    raise ReadError(
                        f"Tecplot: zone {idx + 1} connectivity is truncated"
                    )
                rows.append([int(t) for t in _tokens(self.lines[li])[:npc]])
                li += 1
            conn = np.array(rows, dtype=np.int64).reshape(-1, npc) - 1
        return cols, conn


def load(filename):
    """Parses an ASCII file: (variables, zones, source)."""
    with open_file(filename, "r") as f:
        raw = f.read().splitlines()
    lines = [s for s in (r.strip() for r in raw) if s and not s.startswith("#")]

    variables = []
    zones = []
    i = 0
    while i < len(lines):
        u = lines[i].upper()
        if u.startswith("VARIABLES"):
            joined = lines[i]
            while i + 1 < len(lines) and lines[i + 1].startswith('"'):
                i += 1
                joined += " " + lines[i]
            variables = _variables(joined)
            i += 1
            continue
        if not u.startswith("ZONE"):
            i += 1
            continue
        joined = lines[i]
        while i + 1 < len(lines) and not _is_float(_tokens(lines[i + 1])[0]):
            i += 1
            joined += " " + lines[i]
        z = Zone()
        z.data_start = i + 1
        z.varloc = ""
        fields = {}
        tk = header_tokens(joined)
        k = 1
        while k + 2 < len(tk):
            if tk[k + 1] != "=":
                k += 1
                continue
            key = tk[k].upper()
            val = tk[k + 2]
            if key in (
                "NODES",
                "N",
                "ELEMENTS",
                "E",
                "DATAPACKING",
                "ZONETYPE",
                "F",
                "ET",
                "NV",
                "FACES",
                "TOTALNUMFACENODES",
                "NUMCONNECTEDBOUNDARYFACES",
                "TOTALNUMBOUNDARYCONNECTIONS",
            ):
                fields[key] = val
            elif key in ("I", "J", "K"):
                fields[key] = val
            elif key == "T":
                z.title = _unquote(val)
            elif key == "SOLUTIONTIME":
                z.solution_time = float(val)
                z.has_solution_time = True
            elif key == "STRANDID":
                z.strand = int(val)
                z.has_strand = True
            elif key == "VARLOCATION":
                z.varloc = val.replace(" ", "")
            elif key == "VARSHARELIST":
                for entry in split_entries(val[1:-1]):
                    if "=" in entry:
                        rng, src = entry.split("=", 1)
                        src = int(src) - 1
                    else:  # no zone: the previous one
                        rng, src = entry, len(zones) - 1
                    for idx in parse_ranges(rng):
                        z.var_share[idx] = src
            elif key == "PASSIVEVARLIST":
                z.passive.update(parse_ranges(val[1:-1]))
            elif key == "CONNECTIVITYSHAREZONE":
                z.conn_share = int(val) - 1
            k += 3
        if not variables:
            raise ReadError("Tecplot: no VARIABLES")
        _zone_kind(z, fields, len(variables))
        if z.ordered:
            z.ijk = tuple(int(fields.get(key, 1)) for key in ("I", "J", "K"))
            z.finish()
        else:
            try:
                z.num_nodes = int(fields.get("NODES", fields.get("N")))
                z.num_cells = int(fields.get("ELEMENTS", fields.get("E")))
            except (TypeError, ValueError):
                raise ReadError("Tecplot: an FE zone needs NODES and ELEMENTS")
        if z.is_poly:

            def count(key, required, default):
                if key not in fields:
                    if required:
                        raise ReadError(f"Tecplot: a {z.type_name} zone needs {key}")
                    return default
                try:
                    return int(fields[key])
                except ValueError:
                    raise ReadError(f"Tecplot: bad {key}")

            z.num_faces = count("FACES", True, 0)
            z.total_face_nodes = count(
                "TOTALNUMFACENODES", z.is_polyhedron, 2 * z.num_faces
            )
            z.num_boundary_faces = count("NUMCONNECTEDBOUNDARYFACES", False, 0)
            z.num_boundary_conns = count("TOTALNUMBOUNDARYCONNECTIONS", False, 0)
            if not z.block:
                raise ReadError(
                    f"Tecplot: a {z.type_name} zone must be DATAPACKING=BLOCK"
                )
        want = _data_tokens(z, len(variables))
        li, got = z.data_start, 0
        while got < want and li < len(lines):
            got += len(_tokens(lines[li]))
            li += 1
        if not z.ordered and z.conn_share < 0:
            if z.is_poly:  # the face map is a token stream
                want, seen = _face_map_tokens(z), 0
                while seen < want and li < len(lines):
                    seen += len(_tokens(lines[li]))
                    li += 1
            else:
                li += z.num_cells
        zones.append(z)
        i = li
    if not variables:
        raise ReadError("Tecplot: no VARIABLES")
    if not zones:
        raise ReadError("Tecplot: no ZONE")
    return variables, zones, AsciiSource(lines, zones, len(variables))
