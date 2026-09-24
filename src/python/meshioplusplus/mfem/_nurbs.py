"""
MFEM NURBS meshes (``MFEM NURBS mesh v1.0``/``v1.1``, twin of ``mf_nurbs_*``
in ``mfem.cpp``): the patch topology, its knot vectors and control points,
numbered the way MFEM's ``NURBSExtension`` numbers them (``mesh/nurbs.cpp``,
BSD-3-Clause), so the knot-span elements, their vertices and the control
points of every patch are MFEM's.

The global form keeps every control point once, in MFEM's order: the
topological vertices, then the interior points of every edge, face and patch.
Each patch reaches its own through ``NURBSPatchMap``: the patch's vertices,
its edges and faces with their orientations, and its interior block. The
``patches`` form lists every patch's control points itself (homogeneous
``x w, y w, w`` or ``controlpoints_cartesian``), merged into the same
numbering so grid functions index them alike.
"""

import numpy as np

from .._exceptions import ReadError

# MFEM's local edges and faces (fem/geom.cpp)
QUAD_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0)]
HEX_EDGES = [
    (0, 1),
    (1, 2),
    (3, 2),
    (0, 3),
    (4, 5),
    (5, 6),
    (7, 6),
    (4, 7),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7),
]
HEX_FACES = [
    (3, 2, 1, 0),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
    (4, 5, 6, 7),
]
# geometry of a patch / of its knot-span elements by dimension
GEOM_OF_DIM = {1: 1, 2: 3, 3: 5}
BDR_GEOM_OF_DIM = {1: 0, 2: 1, 3: 3}


def _flip_sign(k):
    return -1 - k


def _unsign(k):
    return k if k >= 0 else -1 - k


class KnotVector:
    def __init__(self, order, knots):
        self.order = order
        self.knots = np.asarray(knots, dtype=np.float64)
        self.ncp = len(self.knots) - order - 1
        # the knot spans of nonzero length: the elements
        self.spans = [
            s for s in range(order, self.ncp) if self.knots[s + 1] > self.knots[s]
        ]

    @property
    def ne(self):
        return len(self.spans)

    def flipped(self):
        k = self.knots
        return KnotVector(self.order, (k[0] + k[-1]) - k[::-1])

    def basis(self, span, r):
        """The ``order + 1`` B-spline values at reference ``r`` in [0, 1] of
        knot span ``span`` (Cox-de Boor), for control points span-order..span."""
        k, p = self.knots, self.order
        u = k[span] + r * (k[span + 1] - k[span])
        n = [1.0] + [0.0] * p
        left = [0.0] * (p + 1)
        right = [0.0] * (p + 1)
        for j in range(1, p + 1):
            left[j] = u - k[span + 1 - j]
            right[j] = k[span + j] - u
            saved = 0.0
            for q in range(j):
                t = n[q] / (right[q + 1] + left[j - q])
                n[q] = saved + right[q + 1] * t
                saved = left[j - q] * t
            n[j] = saved
        return n


def _read_knot(lex):
    line = lex.line()
    order = lex.int("a knot vector order")
    ncp = lex.int("a control point count")
    if order < 0 or ncp < order + 1:
        lex.fail(f"knot vector of order {order} with {ncp} control points", line)
    knots = [lex.real("a knot") for _ in range(ncp + order + 1)]
    if any(b < a for a, b in zip(knots, knots[1:])):
        lex.fail("knots out of order", line)
    return KnotVector(order, knots)


def _quad_orientation(base, test):
    """``Mesh::GetQuadOrientation``."""
    i = next((i for i in range(4) if test[i] == base[0]), 3)
    return 2 * i if test[(i + 1) % 4] == base[1] else 2 * i + 1


def _or1d(n, big_n, o):
    return n if o > 0 else big_n - 1 - n


def _or2d(n1, n2, big_n1, big_n2, o):
    return [
        n1 + n2 * big_n1,
        n2 + n1 * big_n2,
        n2 + (big_n1 - 1 - n1) * big_n2,
        (big_n1 - 1 - n1) + n2 * big_n1,
        (big_n1 - 1 - n1) + (big_n2 - 1 - n2) * big_n1,
        (big_n2 - 1 - n2) + (big_n1 - 1 - n1) * big_n2,
        (big_n2 - 1 - n2) + n1 * big_n2,
        n1 + (big_n2 - 1 - n2) * big_n1,
    ][o]


def _f(n, big_n):
    return 0 if n < 0 else (2 if n >= big_n else 1)


class _PatchMap:
    """``NURBSPatchMap``: patch lattice index -> global number, for the mesh
    vertices (``space=False``) or the control points (``space=True``)."""

    def __init__(self, verts, edges, oedge, faces, oface, poff, counts, opatch):
        self.verts, self.edges, self.oedge = verts, edges, oedge
        self.faces, self.oface = faces, oface
        self.poff, self.opatch = poff, opatch
        self.n = counts  # interior counts per direction

    def ec(self, e, n, big_n, s=1):
        return self.edges[e] + _or1d(n, big_n, s * self.oedge[e])

    def fc(self, f, m, n, big_m, big_n):
        return self.faces[f] + _or2d(m, n, big_m, big_n, self.oface[f])

    def __call__(self, *ijk):
        if len(ijk) == 1:
            i1 = ijk[0] - 1
            c = _f(i1, self.n[0])
            if c == 0:
                return self.verts[0]
            if c == 2:
                return self.verts[1]
            return self.poff + _or1d(i1, self.n[0], self.opatch)
        if len(ijk) == 2:
            i1, j1 = ijk[0] - 1, ijk[1] - 1
            big_i, big_j = self.n[0], self.n[1]
            c = 3 * _f(j1, big_j) + _f(i1, big_i)
            return [
                lambda: self.verts[0],
                lambda: self.ec(0, i1, big_i),
                lambda: self.verts[1],
                lambda: self.ec(3, j1, big_j, -1),
                lambda: self.poff + _or2d(i1, j1, big_i, big_j, self.opatch),
                lambda: self.ec(1, j1, big_j),
                lambda: self.verts[3],
                lambda: self.ec(2, i1, big_i, -1),
                lambda: self.verts[2],
            ][c]()
        i1, j1, k1 = ijk[0] - 1, ijk[1] - 1, ijk[2] - 1
        big_i, big_j, big_k = self.n
        c = 3 * (3 * _f(k1, big_k) + _f(j1, big_j)) + _f(i1, big_i)
        return [
            lambda: self.verts[0],
            lambda: self.ec(0, i1, big_i),
            lambda: self.verts[1],
            lambda: self.ec(3, j1, big_j),
            lambda: self.fc(0, i1, big_j - 1 - j1, big_i, big_j),
            lambda: self.ec(1, j1, big_j),
            lambda: self.verts[3],
            lambda: self.ec(2, i1, big_i),
            lambda: self.verts[2],
            lambda: self.ec(8, k1, big_k),
            lambda: self.fc(1, i1, k1, big_i, big_k),
            lambda: self.ec(9, k1, big_k),
            lambda: self.fc(4, big_j - 1 - j1, k1, big_j, big_k),
            lambda: self.poff + big_i * (big_j * k1 + j1) + i1,
            lambda: self.fc(2, j1, k1, big_j, big_k),
            lambda: self.ec(11, k1, big_k),
            lambda: self.fc(3, big_i - 1 - i1, k1, big_i, big_k),
            lambda: self.ec(10, k1, big_k),
            lambda: self.verts[4],
            lambda: self.ec(4, i1, big_i),
            lambda: self.verts[5],
            lambda: self.ec(7, j1, big_j),
            lambda: self.fc(5, i1, j1, big_i, big_j),
            lambda: self.ec(5, j1, big_j),
            lambda: self.verts[7],
            lambda: self.ec(6, i1, big_i),
            lambda: self.verts[6],
        ][c]()


class Nurbs:
    """A parsed NURBS mesh: topology, knot vectors, control points."""

    def __init__(self, lex, filename, header):
        self.filename = filename
        self.dim = -1
        self.patches = []  # (attr, geom, verts, line)
        self.bpatches = []
        self.edge_rows = []  # (kv, v0, v1) as the file gives them
        self.nv = 0
        self.kvs = None
        self.patch_data = None  # the patches form: [(kvs, sdim, homogeneous cps)]
        self.weights = None
        self.space = None
        self.nodes = None
        self._parse(lex, header)
        self._topology()
        self._offsets()
        self._parse_tail(lex)

    # --- parsing -----------------------------------------------------------------
    def _parse(self, lex, header):
        from ._mfem import _read_elements

        def section(name):
            text, line, _ = lex.next(f"'{name}'")
            if text != name:
                lex.fail(f"expected '{name}', found '{text}'", line)
            return line

        section("dimension")
        self.dim = lex.int("a dimension")
        if self.dim < 1 or self.dim > 3:
            lex.fail(f"dimension {self.dim} (1, 2 or 3)", lex.line())
        section("elements")
        self.patches = _read_elements(lex, "patch")
        section("boundary")
        self.bpatches = _read_elements(lex, "boundary patch")
        for attr, geom, verts, line in self.patches:
            if geom != GEOM_OF_DIM[self.dim]:
                lex.fail(f"a patch of geometry {geom} in a {self.dim}-D mesh", line)
        for attr, geom, verts, line in self.bpatches:
            if geom != BDR_GEOM_OF_DIM[self.dim]:
                lex.fail(f"a boundary patch of geometry {geom}", line)
        section("edges")
        for _ in range(lex.int("an edge count")):
            self.edge_rows.append(
                (lex.int("a knot vector"), lex.int("a vertex"), lex.int("a vertex"))
            )
        section("vertices")
        self.nv = lex.int("a vertex count")
        for group in (self.patches, self.bpatches):
            for _, _, verts, line in group:
                if any(v < 0 or v >= self.nv for v in verts):
                    lex.fail(f"vertex out of range ({self.nv} vertices)", line)
        if not self.edge_rows and self.dim > 1:
            raise ReadError(
                f"MFEM NURBS mesh: {self.filename} has no edges section (edges "
                "MFEM derives itself are not supported)"
            )
        text, line, _ = lex.next("'knotvectors' or 'patches'")
        if text == "knotvectors":
            self.kvs = [_read_knot(lex) for _ in range(lex.int("a count"))]
            if header.endswith("v1.1"):
                self._skip_spacing(lex)
        elif text == "patches":
            self.patch_data = []
            for _ in self.patches:
                section("knotvectors")
                kvs = [_read_knot(lex) for _ in range(lex.int("a count"))]
                if len(kvs) != self.dim:
                    lex.fail(f"a patch with {len(kvs)} knot vectors", line)
                section("dimension")
                sdim = lex.int("a space dimension")
                t, tline, _ = lex.next("control points")
                if t not in (
                    "controlpoints",
                    "controlpoints_homogeneous",
                    "controlpoints_cartesian",
                ):
                    lex.fail(f"expected control points, found '{t}'", tline)
                n = int(np.prod([kv.ncp for kv in kvs]))
                cps = np.array(
                    [lex.real("a coordinate") for _ in range(n * (sdim + 1))]
                ).reshape(n, sdim + 1)
                if t == "controlpoints_cartesian":
                    cps[:, :sdim] *= cps[:, sdim : sdim + 1]
                self.patch_data.append((kvs, sdim, cps))
        else:
            lex.fail(f"expected 'knotvectors' or 'patches', found '{text}'", line)

    def _parse_tail(self, lex):
        """What follows the knot vectors or patches, once the control points
        are counted: the weights and the nodes."""
        while not lex.at_end():
            text, line, _ = lex.next("a section")
            if text in ("mesh_elements", "periodic"):
                raise ReadError(
                    f"MFEM NURBS mesh: the '{text}' section of {self.filename} is not "
                    "supported"
                )
            if text == "weights":
                if self.patch_data is not None:
                    lex.fail("weights in a mesh whose patches carry their own", line)
                self.weights = np.array(
                    [lex.real("a weight") for _ in range(self.num_dofs)]
                )
            elif text in ("unitweights", "autoweights"):
                self.weights = np.ones(self.num_dofs)
            elif text == "FiniteElementSpace":
                self._nodes(lex, line, False)
            elif text == "MFEM" and lex.peek() == "FiniteElementSpace":
                lex.next("FiniteElementSpace")
                version, vline, _ = lex.next("a version")
                if version != "v1.0":
                    lex.fail(f"FiniteElementSpace version '{version}'", vline)
                self._nodes(lex, line, True)
            elif text == "mfem_mesh_end":
                return
            else:
                lex.fail(f"unexpected '{text}'", line)

    def _nodes(self, lex, line, versioned):
        from ._mfem import _read_space_body, _read_space_end

        space = _read_space_body(lex, line)
        if space is None or space.kind != "nurbs":
            lex.fail("NURBS mesh nodes outside a NURBS space", line)
        if versioned:
            _read_space_end(lex)
        self.space = space
        self.nodes = np.asarray(lex.reals(), dtype=np.float64)

    @staticmethod
    def _skip_spacing(lex):
        text = lex.peek()
        if text == "refinements":
            lex.next("refinements")
            lex.reals()
            text = lex.peek()
        if text == "knotvector_refinements":
            lex.next("knotvector_refinements")
            lex.reals()
            text = lex.peek()
        if text != "spacing":
            return
        lex.next("spacing")
        for _ in range(lex.int("a spacing count")):
            lex.int("a knot vector")
            lex.int("a spacing type")
            ni = lex.int("a parameter count")
            nr = lex.int("a parameter count")
            for _ in range(ni):
                lex.int("a parameter")
            for _ in range(nr):
                lex.real("a parameter")

    # --- topology (patchTopo) ----------------------------------------------------
    def _topology(self):
        dim = self.dim
        np_ = len(self.patches)
        if dim == 1:
            if len(self.edge_rows) != np_:
                raise ReadError(
                    f"MFEM NURBS mesh: {len(self.edge_rows)} edges for {np_} patches"
                )
            self.edge_to_ukv = [
                kv if a <= b else _flip_sign(kv) for kv, a, b in self.edge_rows
            ]
        else:
            self.edge_of = {}
            self.edge_to_ukv = []
            for e, (kv, a, b) in enumerate(self.edge_rows):
                self.edge_of[(min(a, b), max(a, b))] = e
                self.edge_to_ukv.append(kv if a <= b else _flip_sign(kv))
        self.pedges = []  # per patch: (edges, orientations)
        for _, _, v, line in self.patches:
            if dim == 1:
                self.pedges.append(None)
                continue
            local = QUAD_EDGES if dim == 2 else HEX_EDGES
            es, os_ = [], []
            for a, b in local:
                key = (min(v[a], v[b]), max(v[a], v[b]))
                if key not in self.edge_of:
                    raise ReadError(
                        f"MFEM NURBS mesh: the edge {v[a]}-{v[b]} of the patch on line "
                        f"{line} is not in the edges section"
                    )
                es.append(self.edge_of[key])
                os_.append(1 if v[a] < v[b] else -1)
            self.pedges.append((es, os_))
        self.pfaces = [None] * np_
        self.face_verts = []
        if dim == 3:
            index = {}
            for p, (_, _, v, _) in enumerate(self.patches):
                fs, os_ = [], []
                for fv in HEX_FACES:
                    verts = tuple(v[c] for c in fv)
                    key = tuple(sorted(verts))
                    if key not in index:
                        index[key] = len(self.face_verts)
                        self.face_verts.append(verts)
                        fs.append(index[key])
                        os_.append(0)
                    else:
                        f = index[key]
                        fs.append(f)
                        os_.append(_quad_orientation(self.face_verts[f], verts))
                self.pfaces[p] = (fs, os_)
            self.face_index = index
        self._boundary_patches()
        self._knot_vectors()

    def _boundary_patches(self):
        """``FinalizeTopology`` and ``CheckBdrElementOrientation``: without a
        boundary section, the faces of one patch become the boundary
        (attribute 1); a boundary patch is turned to run as the face does in
        the first patch holding it."""
        dim = self.dim
        # every face (a vertex, an edge, a quadrilateral): its vertices as the
        # first patch holding it has them, and how many patches hold it
        stored, count = {}, {}
        order = []
        for _, _, v, _ in self.patches:
            if dim == 1:
                faces = [(v[0],), (v[1],)]
            elif dim == 2:
                faces = [(v[a], v[b]) for a, b in QUAD_EDGES]
            else:
                faces = [tuple(v[c] for c in fv) for fv in HEX_FACES]
            for fv in faces:
                key = tuple(sorted(fv))
                if key not in stored:
                    stored[key] = fv
                    count[key] = 0
                    order.append(key)
                count[key] += 1
        if not self.bpatches:
            if dim == 1:
                order.sort()  # 1-D faces are the vertices, numbered as they are
            elif dim == 2:
                order.sort(key=lambda k: self.edge_of[k])
            geom = BDR_GEOM_OF_DIM[dim]
            self.bpatches = [
                (1, geom, list(stored[k]), 0) for k in order if count[k] == 1
            ]
            return
        fixed = []
        for attr, geom, verts, line in self.bpatches:
            key = tuple(sorted(verts))
            fv = stored.get(key)
            verts = list(verts)
            if fv is not None and count[key] == 1:
                if dim == 2 and verts[0] != fv[0]:
                    verts[0], verts[1] = verts[1], verts[0]
                elif dim == 3 and _quad_orientation(fv, verts) % 2:
                    verts[0], verts[2] = verts[2], verts[0]
            fixed.append((attr, geom, verts, line))
        self.bpatches = fixed

    def _kv_sign(self, e):
        return 1 if self.edge_to_ukv[e] >= 0 else -1

    def _kv_of_edge(self, e):
        return self.kvs[_unsign(self.edge_to_ukv[e])]

    def _dir_edges(self, p):
        if self.dim == 1:
            return [p]
        es = self.pedges[p][0]
        return [es[0], es[1]] if self.dim == 2 else [es[0], es[3], es[8]]

    def _kv_dir(self, p):
        """``NURBSExtension::CheckKVDirection``."""
        if self.dim == 1:
            return [self._kv_sign(p)]
        pv = self.patches[p][2]
        pairs = [(pv[0], pv[1]), (pv[0], pv[3])]
        if self.dim == 3:
            pairs.append((pv[0], pv[4]))
        kvdir = [0] * self.dim
        for e in self.pedges[p][0]:
            _, a, b = self.edge_rows[e]
            ks = self._kv_sign(e)
            for d, (x, y) in enumerate(pairs):
                if (a, b) == (x, y):
                    kvdir[d] = ks
                elif (a, b) == (y, x):
                    kvdir[d] = -ks
        return kvdir

    def _knot_vectors(self):
        np_ = len(self.patches)
        if self.patch_data is not None:
            nkv = 1 + max(_unsign(k) for k in self.edge_to_ukv)
            self.kvs = [None] * nkv
            for p in range(np_):
                kvdir = self._kv_dir(p)
                for d, e in enumerate(self._dir_edges(p)):
                    k = _unsign(self.edge_to_ukv[e])
                    if self.kvs[k] is None:
                        kv = self.patch_data[p][0][d]
                        self.kvs[k] = kv.flipped() if kvdir[d] == -1 else kv
            if any(kv is None for kv in self.kvs):
                raise ReadError("MFEM NURBS mesh: a knot vector no patch defines")
        for k in self.edge_to_ukv:
            if _unsign(k) >= len(self.kvs):
                raise ReadError(
                    f"MFEM NURBS mesh: knot vector {_unsign(k)} is not defined"
                )
        self.compr = []
        for p in range(np_):
            kvdir = self._kv_dir(p)
            row = []
            for d, e in enumerate(self._dir_edges(p)):
                kv = self._kv_of_edge(e)
                row.append(kv.flipped() if kvdir[d] == -1 else kv)
            self.compr.append(row)

    # --- numbering (GenerateOffsets) ---------------------------------------------
    def _offsets(self):
        dim = self.dim
        for kind in ("mesh", "space"):
            count = (
                (lambda kv: kv.ne - 1) if kind == "mesh" else (lambda kv: kv.ncp - 2)
            )
            n = self.nv
            eoff = []
            if dim > 1:
                for e in range(len(self.edge_rows)):
                    eoff.append(n)
                    n += count(self._kv_of_edge(e))
            foff = []
            for fv in self.face_verts:
                foff.append(n)
                e0 = self.edge_of[(min(fv[0], fv[1]), max(fv[0], fv[1]))]
                e1 = self.edge_of[(min(fv[1], fv[2]), max(fv[1], fv[2]))]
                n += count(self._kv_of_edge(e0)) * count(self._kv_of_edge(e1))
            poff = []
            for p in range(len(self.patches)):
                poff.append(n)
                size = 1
                for e in self._dir_edges(p):
                    size *= count(self._kv_of_edge(e))
                n += size
            setattr(self, f"{kind}_offsets", (eoff, foff, poff, n))
        self.num_vertices = self.mesh_offsets[3]
        self.num_dofs = self.space_offsets[3]

    def patch_map(self, p, space):
        eoff, foff, poff, _ = self.space_offsets if space else self.mesh_offsets
        kvs = self.compr[p]
        counts = [(kv.ncp - 2) if space else (kv.ne - 1) for kv in kvs]
        verts = list(self.patches[p][2])
        edges, oedge, faces, oface = [], [], [], []
        if self.dim > 1:
            edges = [eoff[e] for e in self.pedges[p][0]]
            oedge = self.pedges[p][1]
        if self.dim == 3:
            faces = [foff[f] for f in self.pfaces[p][0]]
            oface = self.pfaces[p][1]
        return _PatchMap(verts, edges, oedge, faces, oface, poff[p], counts, 0)

    def bdr_patch_map(self, b):
        """``SetBdrPatchVertexMap``: the map and the knot-vector orientations."""
        eoff, foff, _, _ = self.mesh_offsets
        verts = list(self.bpatches[b][2])
        dim = self.dim
        if dim == 1:
            return _PatchMap(verts, [], [], [], [], 0, [0], 0), [1], []
        if dim == 2:
            a, c = verts
            e = self.edge_of.get((min(a, c), max(a, c)))
            if e is None:
                raise ReadError(
                    f"MFEM NURBS mesh: boundary patch {a}-{c} is not an edge of the mesh"
                )
            o = 1 if a < c else -1
            kv = self._kv_of_edge(e)
            okv = o if self.edge_to_ukv[e] >= 0 else -o
            m = _PatchMap(verts, [], [], [], [], eoff[e], [kv.ne - 1], o)
            return m, [okv], [kv]
        es, os_, kvs, okvs = [], [], [], []
        for x, y in QUAD_EDGES:
            a, c = verts[x], verts[y]
            e = self.edge_of.get((min(a, c), max(a, c)))
            if e is None:
                raise ReadError(
                    "MFEM NURBS mesh: a boundary patch edge is not an edge of the mesh"
                )
            es.append(e)
            os_.append(1 if a < c else -1)
        for d in range(2):
            kvs.append(self._kv_of_edge(es[d]))
            okvs.append(os_[d] if self.edge_to_ukv[es[d]] >= 0 else -os_[d])
        f = self.face_index.get(tuple(sorted(verts)))
        if f is None:
            raise ReadError(
                "MFEM NURBS mesh: a boundary patch is not a face of the mesh"
            )
        opatch = _quad_orientation(self.face_verts[f], verts)
        m = _PatchMap(
            verts,
            [eoff[e] for e in es],
            os_,
            [],
            [],
            foff[f],
            [kvs[0].ne - 1, kvs[1].ne - 1],
            opatch,
        )
        return m, okvs, kvs

    # --- the knot-span mesh -------------------------------------------------------
    def elements(self):
        """(attr, geom, verts, line, patch, spans) per knot-span element, in
        MFEM's order."""
        out = []
        geom = GEOM_OF_DIM[self.dim]
        for p, (attr, _, _, line) in enumerate(self.patches):
            m = self.patch_map(p, False)
            kvs = self.compr[p]
            if self.dim == 1:
                for i in range(kvs[0].ne):
                    out.append((attr, geom, [m(i), m(i + 1)], line, p, (i,)))
            elif self.dim == 2:
                for j in range(kvs[1].ne):
                    for i in range(kvs[0].ne):
                        v = [m(i, j), m(i + 1, j), m(i + 1, j + 1), m(i, j + 1)]
                        out.append((attr, geom, v, line, p, (i, j)))
            else:
                for k in range(kvs[2].ne):
                    for j in range(kvs[1].ne):
                        for i in range(kvs[0].ne):
                            v = [
                                m(i, j, k),
                                m(i + 1, j, k),
                                m(i + 1, j + 1, k),
                                m(i, j + 1, k),
                                m(i, j, k + 1),
                                m(i + 1, j, k + 1),
                                m(i + 1, j + 1, k + 1),
                                m(i, j + 1, k + 1),
                            ]
                            out.append((attr, geom, v, line, p, (i, j, k)))
        return out

    def boundary(self):
        out = []
        geom = BDR_GEOM_OF_DIM[self.dim]
        for b, (attr, _, _, line) in enumerate(self.bpatches):
            m, okv, kvs = self.bdr_patch_map(b)
            if self.dim == 1:
                out.append((attr, geom, [m(0)], line))
            elif self.dim == 2:
                nx = kvs[0].ne
                for i in range(nx):
                    i_ = i if okv[0] >= 0 else nx - 1 - i
                    out.append((attr, geom, [m(i_), m(i_ + 1)], line))
            else:
                nx, ny = kvs[0].ne, kvs[1].ne
                for j in range(ny):
                    j_ = j if okv[1] >= 0 else ny - 1 - j
                    for i in range(nx):
                        i_ = i if okv[0] >= 0 else nx - 1 - i
                        v = [m(i_, j_), m(i_ + 1, j_), m(i_ + 1, j_ + 1), m(i_, j_ + 1)]
                        out.append((attr, geom, v, line))
        return out

    # --- control points -----------------------------------------------------------
    def control_points(self):
        """The weights and Cartesian control points (num_dofs x sdim) in MFEM's
        global numbering."""
        if self.patch_data is not None:
            sdim = self.patch_data[0][1]
            w = np.full(self.num_dofs, np.nan)
            x = np.full((self.num_dofs, sdim), np.nan)
            for p, (kvs, psdim, cps) in enumerate(self.patch_data):
                if psdim != sdim:
                    raise ReadError("MFEM NURBS mesh: patches of different dimensions")
                if [kv.ncp for kv in kvs] != [kv.ncp for kv in self.compr[p]]:
                    raise ReadError(
                        "MFEM NURBS mesh: a patch disagrees with its knot vectors"
                    )
                m = self.patch_map(p, True)
                for t, ijk in enumerate(_tensor([kv.ncp for kv in kvs])):
                    g = m(*ijk)
                    w[g] = cps[t, sdim]
                    x[g] = cps[t, :sdim] / cps[t, sdim]
            return w, x
        if self.space is None:
            raise ReadError(f"MFEM NURBS mesh: {self.filename} has no control points")
        from ._mfem import _field_table

        vdim = self.space.vdim
        if len(self.nodes) != self.num_dofs * vdim:
            raise ReadError(
                f"MFEM NURBS mesh: {len(self.nodes)} node values for "
                f"{self.num_dofs} control points of dimension {vdim}"
            )
        x = _field_table(self.nodes, self.space, vdim)
        w = self.weights if self.weights is not None else np.ones(self.num_dofs)
        return w, x

    def evaluate(self, weights, table, element, refs):
        """``table`` (num_dofs x ncomp) at the reference points ``refs`` of a
        knot-span element: the rational (NURBS) interpolant."""
        _, _, _, _, p, spans = element
        kvs = self.compr[p]
        m = self.patch_map(p, True)
        out = np.zeros((len(refs), table.shape[1]))
        for t, r in enumerate(refs):
            per_dir = []
            for d, kv in enumerate(kvs):
                s = kv.spans[spans[d]]
                per_dir.append((s - kv.order, kv.basis(s, r[d])))
            num = np.zeros(table.shape[1])
            den = 0.0
            for idx in _tensor([kv.order + 1 for kv in kvs]):
                n = 1.0
                cp = []
                for d, a in enumerate(idx):
                    n *= per_dir[d][1][a]
                    cp.append(per_dir[d][0] + a)
                g = m(*cp)
                nw = n * weights[g]
                num += nw * table[g]
                den += nw
            out[t] = num / den
        return out


def _tensor(counts):
    """Lexicographic multi-indices, the first fastest."""
    if len(counts) == 1:
        return [(i,) for i in range(counts[0])]
    if len(counts) == 2:
        return [(i, j) for j in range(counts[1]) for i in range(counts[0])]
    return [
        (i, j, k)
        for k in range(counts[2])
        for j in range(counts[1])
        for i in range(counts[0])
    ]
