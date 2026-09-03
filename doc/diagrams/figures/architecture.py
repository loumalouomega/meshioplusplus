"""Architecture and data-model figures."""

from diaglib import palette as P
from diaglib.layout import spread
from diaglib.svg import Canvas
from figures.cells import tables


def architecture():
    c = Canvas(
        960,
        640,
        "meshio++ architecture: one C++ core, six language surfaces, and the tools built on them",
    )
    # consumers
    consumers = [
        (["Python CLI", "meshioplusplus"], P.PYTHON, 0),
        (["native CLI", "(no Python)"], P.CORE, 3),
        (["MCP server", "57 tools"], P.PYTHON, 0),
        (["browser viewer +", "dataset manager"], P.WASM, 2),
        (["Polyscope", "viewer"], P.PYTHON, 0),
        (["Blender add-on", "ParaView plugin"], P.PYTHON, 0),
        (["interop · GPU", "ML · PhysicsNeMo"], P.PYTHON, 0),
        (["Fortran · Julia", "· R programs"], P.CABI, 1),
    ]
    xs = spread(24, 936, len(consumers), 104)
    c.text(
        24,
        36,
        "consumers",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
        fill=P.INK_2,
    )
    for x, (label, colour, _) in zip(xs, consumers):
        c.box(x, 44, 104, 50, label, color=colour, size=P.SIZE_SMALL, weight="600")
    # surfaces
    c.text(
        24,
        150,
        "language surfaces",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
        fill=P.INK_2,
    )
    sx = [24, 282, 552, 762]
    sw = [246, 258, 196, 174]
    c.box(
        sx[0],
        158,
        sw[0],
        78,
        ["Python", "_core (pybind11) + pure-Python fallback"],
        color=P.PYTHON,
        sub="zero-copy numpy at the I/O boundary",
        size=P.SIZE_LABEL,
    )
    c.box(
        sx[1],
        158,
        sw[1],
        78,
        ["C API — libmeshioplusplus", "pure C99 header, SOVERSION 0"],
        color=P.CABI,
        sub="Fortran module · Julia · R ride on it",
        size=P.SIZE_LABEL,
    )
    c.box(
        sx[2],
        158,
        sw[2],
        78,
        ["WebAssembly", "@meshioplusplus/wasm"],
        color=P.WASM,
        sub="embind, MEMFS, threaded + seq",
        size=P.SIZE_LABEL,
    )
    c.box(
        sx[3],
        158,
        sw[3],
        78,
        ["C++ API", "find_package(meshioplusplus)"],
        color=P.CORE,
        sub="or the single-header amalgamation",
        size=P.SIZE_LABEL,
    )
    # consumer -> surface arrows
    for x, (label, colour, target) in zip(xs, consumers):
        tx = sx[target] + sw[target] / 2
        c.elbow(x + 52, 94, tx, 156, stroke=colour, sw=1.4, via="mid")
    # the typed pipeline band
    c.box(
        24,
        262,
        912,
        36,
        "settings pipeline · sequence driver · provenance scope — one typed layer (run_pipeline_steps) that every surface drives, so a settings.json and the viewer's op chain cannot drift",
        color=P.FORMATS,
        size=P.SIZE_SMALL,
        weight="600",
    )
    for x, w in zip(sx, sw):
        c.arrow(x + w / 2, 236, x + w / 2, 260, stroke=P.INK_2, sw=1.4)
    c.arrow(480, 298, 480, 318, stroke=P.INK_2, sw=1.8)
    # the core
    c.rrect(
        24, 320, 912, 288, fill=P.CORE, fill_opacity=0.05, stroke=P.CORE, rx=8, sw=1.8
    )
    c.text(
        40,
        344,
        "one C++ core — meshioplusplus_core_obj (C++20; every optional dependency compiles out when absent)",
        size=P.SIZE_LABEL,
        anchor="start",
        weight="700",
        fill=P.CORE,
    )
    c.box(
        44,
        358,
        280,
        66,
        ["format registry: 43 readable,", "46 writable formats"],
        color=P.FORMATS,
        sub="registry.cpp: name → reader/writer, extension → default",
        size=P.SIZE_SMALL,
    )
    c.box(
        344,
        358,
        280,
        66,
        ["operations layer: 34 mesh", "+ 5 data operations"],
        color=P.CORE,
        sub="uniform-mesh-API only; parallel_for hot loops",
        size=P.SIZE_SMALL,
    )
    c.box(
        644,
        358,
        272,
        66,
        ["detail/ kernels"],
        color=P.CORE,
        sub=[
            "cell_faces · region_remap · marching · spatial_hash",
            "polyhedron · refine_templates · provenance",
        ],
        size=P.SIZE_SMALL,
    )
    c.box(
        44,
        438,
        872,
        26,
        "uniform mesh API (mesh_api.hpp): AssignPoints · AddCellBlock · AddRegion · Points() · CellRange() · PointData(name) …",
        color=P.INK_2,
        size=P.SIZE_SMALL,
        fill_opacity=0.06,
        weight="600",
    )
    c.text(
        40,
        484,
        "mesh backend, exactly one per build (MESHIOPLUSPLUS_MESH_BACKEND):",
        size=P.SIZE_SMALL,
        anchor="start",
        fill=P.INK_2,
    )
    c.box(
        44,
        492,
        280,
        44,
        "MESHIO",
        color=P.PYTHON,
        sub="meshio-mirroring Mesh/CellBlock; the Python wheel",
        size=P.SIZE_SMALL,
    )
    c.box(
        344,
        492,
        280,
        44,
        "NATIVE",
        color=P.WASM,
        sub="canonical Float64/Int64, CSR ragged; the WASM build",
        size=P.SIZE_SMALL,
    )
    c.box(
        644,
        492,
        272,
        44,
        "KRATOS",
        color=P.CABI,
        sub="a Kratos-style ModelPart, materialised lazily",
        size=P.SIZE_SMALL,
    )
    c.box(
        44,
        550,
        872,
        40,
        "optional dependencies: HDF5 · netCDF · zlib / zstd / lz4 · KaHIP · Kokkos · Eigen (submodule) · nlohmann/json (submodule) · gidpost (vendored) · Polyscope (CLI only)",
        color=P.MUTED,
        size=P.SIZE_SMALL,
        fill_opacity=0.08,
        weight="500",
    )
    c.legend(
        24,
        626,
        [
            (P.CORE, "C++"),
            (P.PYTHON, "Python"),
            (P.CABI, "C ABI family"),
            (P.WASM, "WASM / browser"),
            (P.FORMATS, "shared dispatch"),
        ],
    )
    c.label(936, 626, "arrows: who consumes what", anchor="end")
    return c.render()


def mesh_backends():
    c = Canvas(
        900,
        460,
        "Three interchangeable in-memory mesh backends behind one uniform mesh API",
    )
    # left: the callers
    callers = [
        (["format readers / writers", "43 formats, src/cpp/src/formats/"], P.FORMATS),
        (["operations", "34 mesh + 5 data, operations/"], P.CORE),
        (["bindings", "pybind11 · C API · embind"], P.PYTHON),
    ]
    for k, (label, colour) in enumerate(callers):
        y = 70 + k * 100
        c.box(
            24,
            y,
            214,
            70,
            label,
            color=colour,
            size=P.SIZE_LABEL,
            sub_size=P.SIZE_SMALL,
        )
        c.arrow(238, y + 35, 284, y + 35, stroke=colour, sw=1.6)
    c.text(131, 52, "every caller", size=P.SIZE_SMALL, weight="700", fill=P.INK_2)
    # centre: the API
    c.rrect(288, 50, 176, 310, fill=P.INK_2, fill_opacity=0.06, stroke=P.INK_2, rx=8)
    c.text(376, 74, "uniform mesh API", size=P.SIZE_LABEL, weight="700")
    c.text(376, 90, "mesh_api.hpp", size=P.SIZE_SMALL, mono=True, fill=P.INK_2)
    api = [
        "AssignPoints",
        "AddCellBlock",
        "AddPolygonBlock",
        "AddPolyhedronBlock",
        "AddPointData / AddCellData",
        "AddFieldData / AddRegion",
        "AddPropertySet",
        "—",
        "Points() / PointDim()",
        "Cells(i) / CellRange()",
        "PointData(name)",
        "CellData(name, block)",
        "Region(i) / RegionNames()",
        "PointDataNames() (sorted)",
    ]
    for k, s in enumerate(api):
        c.text(
            376,
            114 + k * 17,
            s,
            size=P.SIZE_SMALL,
            mono=s != "—",
            fill=P.INK if s != "—" else P.MUTED,
        )
    c.arrow(464, 205, 508, 205, stroke=P.INK_2, sw=1.8)
    # right: the backends
    backends = [
        (
            "MESHIO (default)",
            P.PYTHON,
            [
                "Mesh / CellBlock mirroring the Python data model; dtype-erased NDArrays",
                "required by the Python extension, which views numpy memory into it",
            ],
        ),
        (
            "NATIVE",
            P.WASM,
            [
                "canonical Float64 points, Int64 connectivity, CSR ragged blocks",
                "lazy whole-mesh CSR (GlobalConnectivity()); the WASM build uses it",
            ],
        ),
        (
            "KRATOS",
            P.CABI,
            [
                "stages canonically, materialises a Kratos-style ModelPart lazily",
                "SubModelParts ⇄ named regions; kratos_bridge.hpp → a real Kratos::ModelPart",
            ],
        ),
    ]
    for k, (name, colour, sub) in enumerate(backends):
        y = 60 + k * 100
        c.box(
            512, y, 364, 84, name, color=colour, sub=sub, size=P.SIZE_LABEL, sub_size=10
        )
    c.text(
        694,
        52,
        "MESHIOPLUSPLUS_MESH_BACKEND — exactly one is compiled in",
        size=P.SIZE_SMALL,
        weight="700",
        fill=P.INK_2,
    )
    c.label(
        450,
        392,
        "a format, an operation or a binding written against the uniform API compiles unchanged under all three backends,",
    )
    c.label(
        450,
        408,
        "which is what CI's cpp-tests matrix (MESHIO×STL, MESHIO×OpenMP, NATIVE×OpenMP, KRATOS×OpenMP) checks on every change;",
    )
    c.label(
        450,
        424,
        "NATIVE and KRATOS canonicalise dtypes at ingest within kind (floats → Float64, ints → Int64, never int → float), and every backend sorts data names,",
    )
    c.label(
        450,
        440,
        "so on-disk field order and region canonicalisation are byte-identical across backends and thread counts",
    )
    return c.render()


def mesh_data_model():
    c = Canvas(
        900,
        430,
        "The Mesh data model: points, cell blocks, and data arrays aligned to them",
    )

    # points table
    def table(x, y, title, rows, colour, w=150, row_h=18, mono=True, header=None):
        c.text(
            x,
            y - 8,
            title,
            size=P.SIZE_SMALL,
            anchor="start",
            weight="700",
            fill=colour,
        )
        c.rrect(
            x,
            y,
            w,
            row_h * len(rows),
            fill=colour,
            fill_opacity=0.08,
            stroke=colour,
            rx=3,
            sw=1.2,
        )
        for k, r in enumerate(rows):
            if k:
                c.line(
                    x,
                    y + k * row_h,
                    x + w,
                    y + k * row_h,
                    stroke=colour,
                    sw=0.6,
                    opacity=0.5,
                )
            c.text(x + 6, y + k * row_h + 13, r, size=10, anchor="start", mono=mono)
        return y + row_h * len(rows)

    pts = [
        "0: 0.0 0.0 0.0",
        "1: 1.0 0.0 0.0",
        "2: 0.0 1.0 0.0",
        "3: 1.0 1.0 0.0",
        "4: 2.0 0.0 0.0",
        "5: 2.0 1.0 0.0",
    ]
    table(30, 70, "points  (6, 3) float64", pts, P.INK, w=140)
    tdata = ["0: 20.0", "1: 21.5", "2: 19.0", "3: 22.0", "4: 24.0", "5: 25.5"]
    table(190, 70, 'point_data["T"]  (6,)', tdata, P.DATA, w=90)
    for k in range(6):
        c.line(
            170,
            70 + k * 18 + 9,
            190,
            70 + k * 18 + 9,
            stroke=P.DATA,
            sw=0.8,
            dash="2 2",
        )
    c.label(155, 200, "one row per point", anchor="middle")
    # cells
    c.text(
        330,
        62,
        "cells: [CellBlock(type, (n, k) ints), …]",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
        fill=P.CORE,
    )
    yb = table(
        330,
        84,
        'CellBlock("triangle", (2, 3))',
        ["0: 0 1 2", "1: 1 3 2"],
        P.CORE,
        w=170,
    )
    table(
        330, yb + 30, 'CellBlock("quad", (1, 4))', ["0: 1 4 5 3"], P.CORE, w=170
    )
    # cell_data aligned to blocks
    c.text(
        540,
        62,
        'cell_data["material"]  — one array per block',
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
        fill=P.DATA,
    )
    table(540, 84, "block 0: (2,)", ["0: 1", "1: 1"], P.DATA, w=90)
    table(540, yb + 30, "block 1: (1,)", ["0: 2"], P.DATA, w=90)
    for y in (84 + 9, 84 + 27, yb + 30 + 9):
        c.line(500, y, 540, y, stroke=P.DATA, sw=0.8, dash="2 2")
    # global block-major ruler
    c.text(
        330,
        212,
        "global cell index (block-major, detail/cell_index.hpp): ",
        size=P.SIZE_SMALL,
        anchor="start",
        fill=P.INK_2,
    )
    for k, (label, colour) in enumerate([("0", P.CORE), ("1", P.CORE), ("2", P.CORE)]):
        c.rrect(
            330 + k * 34,
            220,
            30,
            20,
            fill=colour,
            fill_opacity=0.12,
            stroke=colour,
            rx=3,
            sw=1,
        )
        c.text(345 + k * 34, 234, label, size=P.SIZE_SMALL, mono=True)
    c.label(330, 256, "block 0 → 0, 1   ·   block 1 → 2", anchor="start")
    # field_data and regions
    table(
        660,
        84,
        "field_data",
        ['"time": 0.5', '"exodus:time": …'],
        P.REGIONS,
        w=180,
        mono=True,
    )
    c.text(
        660,
        150,
        "regions  (named groups)",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
        fill=P.REGIONS,
    )
    rows = [
        'Region("inlet", point, [0, 2])',
        'Region("solid", cell, [0, 1, 2])',
        'Region("wall", side, [(2, 0)])',
    ]
    c.rrect(
        660,
        158,
        210,
        54,
        fill=P.REGIONS,
        fill_opacity=0.08,
        stroke=P.REGIONS,
        rx=3,
        sw=1.2,
    )
    for k, r in enumerate(rows):
        c.text(666, 171 + k * 18, r, size=10, anchor="start", mono=True)
    c.label(660, 230, "point_sets / cell_sets are dict views over them", anchor="start")
    c.label(660, 246, "side entries are (global cell, local facet)", anchor="start")
    # explanatory text
    lines = [
        "A Mesh is points plus a list of homogeneous cell blocks. Every array describing cells is a LIST with exactly one entry per block, in block order:",
        "the one-array-per-block invariant every operation preserves (an op writing cell_data emits NumCellBlocks() arrays; an op reading one checks the count).",
        "Point data is one array aligned to points; field data is whole-mesh. Regions address cells by their global block-major index, so a region survives",
        "every operation that carries a cell map (detail/region_remap.hpp owns that remapping). Data-name order is observable — it drives most writers'",
        "on-disk field order — so every backend returns sorted names.",
    ]
    for k, s in enumerate(lines):
        c.label(450, 296 + 17 * k, s)
    c.legend(
        30,
        410,
        [
            (P.INK, "geometry"),
            (P.CORE, "connectivity"),
            (P.DATA, "data arrays"),
            (P.REGIONS, "regions and field data"),
        ],
    )
    return c.render()


def regions():
    c = Canvas(
        900,
        520,
        "Named regions: every format's group concept maps onto one Region model of three kinds",
    )
    sources = [
        ("gmsh", "physical group (per dim, tagged)"),
        ("Exodus", "element block · node set · side set"),
        ("Abaqus", "*NSET · *ELSET · *SURFACE"),
        ("MED", "family + GRO group names"),
        ("UNV", "group (Phase 2)"),
        ("Ansys", "component (Phase 2)"),
        ("OpenFOAM", "boundary patch (Phase 2)"),
        ("Kratos", "SubModelPart"),
    ]
    for k, (fmt, what) in enumerate(sources):
        y = 44 + k * 44
        phase2 = "Phase 2" in what
        c.box(
            24,
            y,
            230,
            34,
            [f"{fmt}: {what.replace(' (Phase 2)', '')}"],
            color=P.FORMATS if not phase2 else P.MUTED,
            size=P.SIZE_SMALL,
            weight="600",
            dash="4 3" if phase2 else None,
        )
        c.arrow(
            254, y + 17, 318, 210, stroke=P.FORMATS if not phase2 else P.MUTED, sw=1.2
        )
    # the Region struct
    c.rrect(
        320,
        130,
        220,
        160,
        fill=P.REGIONS,
        fill_opacity=0.10,
        stroke=P.REGIONS,
        rx=8,
        sw=1.6,
    )
    c.text(430, 154, "Region", size=P.SIZE_TITLE, weight="700", fill=P.REGIONS)
    fields = [
        "name: str",
        "kind: point | cell | side",
        "dim: int  (-1 = unspecified)",
        "tag: int  (-1 = none)",
        "entries: int64 array, sorted,",
        "         de-duplicated (canonical)",
    ]
    for k, s in enumerate(fields):
        c.text(332, 178 + k * 17, s, size=P.SIZE_SMALL, anchor="start", mono=True)
    c.arrow(540, 210, 604, 210, stroke=P.REGIONS, sw=1.6)
    # the three kinds
    kinds = [
        ("point", "entries are point indices, shape (n,)", P.INK),
        ("cell", "entries are GLOBAL block-major cell indices, (n,)", P.CORE),
        ("side", "entries are (global cell, local facet) pairs, (n, 2)", P.CABI),
    ]
    for k, (name, what, colour) in enumerate(kinds):
        y = 100 + k * 76
        c.box(
            608,
            y,
            268,
            60,
            [f'kind = "{name}"'],
            color=colour,
            sub=what,
            size=P.SIZE_LABEL,
            sub_size=10,
            mono=False,
        )
    c.label(
        742,
        340,
        "facets numbered per detail/cell_faces.hpp / cell_edges.hpp;",
        anchor="middle",
    )
    c.label(
        742,
        356,
        "side regions have no point_sets / cell_sets equivalent",
        anchor="middle",
    )
    # block-major ruler
    y = 430
    c.text(
        24,
        y - 12,
        "global cell index vs (block, row) — detail/cell_index.hpp, the single owner:",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
        fill=P.INK_2,
    )
    blocks = [
        ("block 0: triangle", 5, P.CORE),
        ("block 1: quad", 3, P.CORE),
        ("block 2: tetra", 4, P.CORE),
    ]
    x = 24
    g = 0
    for label, n, colour in blocks:
        c.text(x, y + 44, label, size=P.SIZE_SMALL, anchor="start", fill=P.INK_2)
        for r in range(n):
            c.rrect(
                x, y, 30, 24, fill=colour, fill_opacity=0.12, stroke=colour, rx=3, sw=1
            )
            c.text(x + 15, y + 16, str(g), size=P.SIZE_SMALL, mono=True)
            c.text(x + 15, y + 34, f"r{r}", size=9, fill=P.MUTED)
            x += 32
            g += 1
        x += 22
    c.label(
        x + 10,
        y + 16,
        "a cell region {3, 7, 9} = (block 0, row 3), (block 1, row 2), (block 2, row 1)",
        anchor="start",
    )
    c.label(
        x + 10,
        y + 34,
        "the same numbering partition_labels and every cell map use",
        anchor="start",
    )
    c.legend(
        24,
        510,
        [
            (P.FORMATS, "maps regions today"),
            (P.MUTED, "Phase 2: not yet mapped"),
            (P.REGIONS, "the one model"),
        ],
    )
    return c.render()


def polyhedra_csr():
    c = Canvas(
        900,
        440,
        "Ragged cell blocks cross every flat binding as CSR arrays: cell_offsets → face_offsets → nodes",
    )
    t = tables()
    # two cells: a tetra as polyhedron4 (nodes 0-3) and a pyramid as polyhedron5 (nodes 4-8)
    cells = []
    for base, offset in (("tetra", 0), ("pyramid", 4)):
        cells.append(
            [tuple(i + offset for i in row[2][: row[1]]) for row in t.cell_faces[base]]
        )
    cell_offsets = [0]
    face_offsets = [0]
    nodes = []
    for faces in cells:
        for f in faces:
            nodes.extend(f)
            face_offsets.append(len(nodes))
        cell_offsets.append(cell_offsets[-1] + len(faces))

    def strip(x, y, title, values, colour, w=22, highlight=(), tag=None):
        c.text(
            x,
            y - 8,
            title,
            size=P.SIZE_SMALL,
            anchor="start",
            weight="700",
            fill=colour,
            mono=True,
        )
        for k, v in enumerate(values):
            hi = k in highlight
            c.rrect(
                x + k * w,
                y,
                w,
                22,
                fill=colour,
                fill_opacity=0.35 if hi else 0.10,
                stroke=colour,
                rx=2,
                sw=1.4 if hi else 0.8,
            )
            c.text(
                x + k * w + w / 2,
                y + 15,
                str(v),
                size=10,
                mono=True,
                weight="700" if hi else None,
            )
            c.text(x + k * w + w / 2, y + 34, str(k), size=8, fill=P.MUTED)
        return x + len(values) * w

    c.caption(
        450,
        34,
        "2-level block (polyhedron): cell 0 = a tetra-shaped polyhedron4 (nodes 0–3), cell 1 = a pyramid-shaped polyhedron5 (nodes 4–8)",
    )
    strip(
        30,
        70,
        "cell_offsets  (num_cells + 1)  — each cell's first face",
        cell_offsets,
        P.CABI,
        w=30,
        highlight=(1, 2),
    )
    strip(
        30,
        140,
        "face_offsets  (num_faces + 1)  — each face's first node",
        face_offsets,
        P.PYTHON,
        w=30,
        highlight=(4, 5),
    )
    strip(
        30,
        210,
        "nodes  (num_nodes)  — every face's node ids, concatenated",
        nodes,
        P.CORE,
        w=22,
        highlight=tuple(range(face_offsets[4], face_offsets[5])),
    )
    # arrows: cell 1, face 0
    c.arrow(30 + 1 * 30 + 15, 92, 30 + 4 * 30 + 15, 138, stroke=P.CABI, sw=1.6)
    c.text(
        140,
        86,
        "cell 1 starts at face 4 (cell_offsets[1])",
        size=P.SIZE_SMALL,
        anchor="start",
        fill=P.CABI,
    )
    c.arrow(
        30 + 4 * 30 + 15,
        162,
        30 + face_offsets[4] * 22 + 11,
        208,
        stroke=P.PYTHON,
        sw=1.6,
    )
    c.text(
        30 + 10 * 30 + 12,
        156,
        f"face 4 = nodes[{face_offsets[4]}..{face_offsets[5]}) — the pyramid's quad base",
        size=P.SIZE_SMALL,
        anchor="start",
        fill=P.PYTHON,
    )
    c.text(
        30,
        270,
        "face f of cell c  =  nodes[ face_offsets[ cell_offsets[c] + f ] .. face_offsets[ cell_offsets[c] + f + 1 ] )",
        size=P.SIZE_SMALL,
        anchor="start",
        mono=True,
    )
    c.label(
        30,
        288,
        "the pyramid's first face is its quad base: 4 node ids, where the tetra's faces have 3 — that is what a rectangular array cannot hold",
        anchor="start",
    )
    # 1-level inset
    c.rrect(
        30,
        310,
        840,
        70,
        fill=P.FORMATS,
        fill_opacity=0.06,
        stroke=P.FORMATS,
        rx=6,
        sw=1.2,
    )
    c.text(
        44,
        330,
        "1-level block (jagged polygon rows): one face per cell, so face_offsets IS the row-offsets array and there is no cell_offsets at all",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="600",
    )
    c.text(
        44,
        348,
        "the C API returns NULL and JS omits the key, rather than synthesising a 0, 1, 2, … identity a caller could mistake for information",
        size=P.SIZE_SMALL,
        anchor="start",
        fill=P.INK_2,
    )
    c.text(
        44,
        366,
        "JS names: data · faceOffsets (rowOffsets when 1-level) · cellOffsets; NATIVE stores the same two arrays under swapped names (mRowOffsets / mFaceOffsets)",
        size=P.SIZE_SMALL,
        anchor="start",
        fill=P.INK_2,
    )
    c.label(
        450,
        410,
        "winding is repaired per cell before any measurement (never required), and a non-planar face is measured by the corner-average fan — see the page text",
    )
    return c.render()


FIGURES = {
    "architecture": architecture,
    "mesh_backends": mesh_backends,
    "mesh_data_model": mesh_data_model,
    "regions": regions,
    "polyhedra_csr": polyhedra_csr,
}
