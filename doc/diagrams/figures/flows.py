"""Flow and decision figures: pipeline, sequences, provenance, interop, data, ML, ABI, roadmap."""

from diaglib import palette as P
from diaglib.layout import spread
from diaglib.svg import Canvas


def pipeline_flow():
    c = Canvas(
        900,
        440,
        "The settings pipeline: one JSON document, validated before reading, run by one typed layer on every surface",
    )
    steps = [
        (["settings.json"], P.FORMATS, "four top-level keys"),
        (["parse", "(strict)"], P.INK_2, "unknown key → error"),
        (["validate", "EVERY step"], P.ERROR, "before any read"),
        (["read", "Input.Options"], P.FORMATS, "→ ReadOptions"),
        (["op₁ … opₙ"], P.CORE, "PipelineStep{Op, …}"),
        (["write", "Output"], P.FORMATS, "registry_write_ex"),
        (["report"], P.DATA, "{steps, warnings}"),
    ]
    xs = spread(24, 876, len(steps), 108)
    for x, (label, colour, note) in zip(xs, steps):
        c.box(x, 110, 108, 56, label, color=colour, size=P.SIZE_LABEL)
        c.label(x + 54, 186, note, size=9.5)
    for a, b in zip(xs, xs[1:]):
        c.arrow(a + 108, 138, b, 138, stroke=P.INK_2, sw=1.6)
    # the sequence branch
    c.box(
        xs[3] - 30,
        30,
        360,
        44,
        ["Input names a multi-step file or a glob, or a sequence key is present?"],
        color=P.CABI,
        size=P.SIZE_SMALL,
        sub="→ the sequence driver runs the same chain once per step, one mesh alive at a time",
        sub_size=9.5,
    )
    c.arrow(xs[3] + 54, 108, xs[3] + 54, 76, stroke=P.CABI, sw=1.4, dash="4 3")
    c.arrow(xs[5] + 54, 76, xs[5] + 54, 108, stroke=P.CABI, sw=1.4, dash="4 3")
    # the single owner
    c.box(
        24,
        230,
        852,
        40,
        "run_pipeline_steps — the single owner of the step dispatch (operations/pipeline.hpp); pipeline_op_table() lists every op and key",
        color=P.CORE,
        size=P.SIZE_LABEL,
        weight="700",
    )
    surfaces = [
        ("both CLIs", "pipeline verb", P.CORE),
        ("Python", "run_pipeline (pure twin)", P.PYTHON),
        ("C / Fortran / Julia / R", "mio_pipeline_run_file/json", P.CABI),
        ("WASM", "runPipeline / convertSurfaceOps", P.WASM),
        ("MCP", "pipeline tool, sandboxed", P.PYTHON),
        ("browser viewer", "its op chain, replayed", P.WASM),
    ]
    sx = spread(24, 876, len(surfaces), 132)
    for x, (name, what, colour) in zip(sx, surfaces):
        c.box(
            x,
            310,
            132,
            52,
            name,
            color=colour,
            sub=what,
            size=P.SIZE_SMALL,
            sub_size=9.5,
        )
        c.arrow(x + 66, 308, x + 66, 272, stroke=colour, sw=1.4)
    c.label(
        450,
        392,
        'vocabulary is PascalCase ({"Op": "ConvertCells", "Mode": "simplexify"}); the viewer\'s camelCase convertSurfaceOps is the same table re-cased mechanically;',
    )
    c.label(
        450,
        408,
        "multi-mesh ops (Merge, Interpolate, ConservativeInterpolate, Split, Diff, UndoGreen) are refused by name, pointing at their CLI verb — v1 is a single-mesh chain;",
    )
    c.label(
        450,
        424,
        "the Python engine is a deliberate pure-Python twin over the public API, pinned mesh-for-mesh against the C++ engine by tests/python/test_pipeline.py",
    )
    return c.render()


def _file(c, x, y, label, colour=P.FORMATS, w=96, h=40, sub=None):
    c.box(
        x,
        y,
        w,
        h,
        label,
        color=colour,
        size=P.SIZE_SMALL,
        sub=sub,
        sub_size=9,
        mono=True,
        rx=4,
    )


def sequences_shapes():
    c = Canvas(
        900,
        480,
        "Sequences: fan-in, fan-out, per-step chains, and the truncation that is refused",
    )
    rows = [
        ("fan-in", "N single-step files → one multi-step file (XDMF)", P.CABI),
        ("fan-out", "one multi-step file → one file per step", P.CABI),
        ("sequence (N → N)", "the operation chain applied to each step", P.CABI),
        ("refused", "a multi-step input aimed at a single-step output", P.ERROR),
    ]
    for k, (name, what, colour) in enumerate(rows):
        y = 40 + k * 100
        c.text(
            24,
            y + 14,
            name,
            size=P.SIZE_LABEL,
            anchor="start",
            weight="700",
            fill=colour,
        )
        c.label(24, y + 30, what, anchor="start")
        c.line(24, y + 82, 876, y + 82, stroke=P.HAIRLINE)
    # row 1: fan-in
    y = 40
    for i in range(3):
        _file(c, 300 + i * 104, y + 4, f"out_000{i}.vtu")
    c.arrow(614, y + 24, 700, y + 24, stroke=P.CABI, sw=1.8)
    _file(c, 706, y + 4, "series.xdmf", colour=P.CABI, w=120, sub="3 steps, one .h5")
    c.label(
        456,
        y + 62,
        "ordering is natural-numeric: out_9 < out_10",
        size=9.5,
        anchor="start",
    )
    # row 2: fan-out
    y = 140
    _file(c, 300, y + 4, "series.xdmf", colour=P.CABI, w=120, sub="3 steps")
    c.arrow(426, y + 24, 512, y + 24, stroke=P.CABI, sw=1.8)
    for i in range(3):
        _file(c, 518 + i * 104, y + 4, f"out_000{i}.vtu")
    c.label(300, y + 62, "Output: out_{step}.vtu", size=9.5, mono=True, anchor="start")
    # row 3: sequence
    y = 240
    for i in range(3):
        _file(c, 300 + i * 104, y + 4, f"in_{i}.vtu")
    c.arrow(614, y + 24, 640, y + 24, stroke=P.CABI, sw=1.8)
    c.box(
        644,
        y + 2,
        110,
        44,
        ["clean → refine", "→ quality"],
        color=P.CORE,
        size=P.SIZE_SMALL,
        sub="per step",
        sub_size=9,
    )
    c.arrow(756, y + 24, 782, y + 24, stroke=P.CABI, sw=1.8)
    for i in range(3):
        _file(c, 786 + i * 0, y + 4 + i * 0, "out_{step}", w=90) if i == 0 else None
    c.label(830, y + 58, "×3, one mesh alive at a time", size=9.5)
    # row 4: refusal
    y = 340
    _file(c, 300, y + 4, "series.xdmf", colour=P.CABI, w=120, sub="3 steps")
    c.arrow(426, y + 24, 512, y + 24, stroke=P.ERROR, sw=1.8, dash="5 3")
    _file(c, 518, y + 4, "out.vtu", colour=P.ERROR)
    c.text(
        640,
        y + 20,
        "✗ error naming the remedy:",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
        fill=P.ERROR,
    )
    c.label(
        640,
        y + 36,
        "use out_{step}.vtu or --time-step N; never step 0 silently",
        anchor="start",
    )
    c.label(
        450,
        448,
        "the mode is inferred from the input's step count and whether the output pattern has {step}/{index}; an explicit Mode only asserts it.",
    )
    c.label(
        450,
        464,
        'time values, in precedence: an explicit list → the file (a series step or field_data["meshio:time"]) → the stem\'s last digit run → the index',
    )
    return c.render()


def provenance_write_paths():
    c = Canvas(
        960,
        540,
        "Where the provenance record lives: a thread-local scope every writer reads, whichever of the write paths reaches it",
    )
    # left column: paths through registry_write_ex
    c.text(
        24,
        40,
        "paths that go through registry_write_ex",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
        fill=P.CORE,
    )
    entries = [
        "native CLI (convert, pipeline)",
        "mio_write_ex (C ABI)",
        "settings pipeline (Output)",
        "sequence driver",
    ]
    for k, e in enumerate(entries):
        c.box(
            24, 52 + k * 40, 236, 32, e, color=P.CORE, size=P.SIZE_SMALL, weight="600"
        )
        c.elbow(260, 68 + k * 40, 300, 240, stroke=P.CORE, sw=1.2, via="h", arrow=False)
    c.box(
        300,
        224,
        180,
        44,
        ["registry_write_ex", "(WriteOptions)"],
        color=P.CORE,
        size=P.SIZE_SMALL,
        mono=True,
    )
    c.arrow(390, 268, 390, 330, stroke=P.CORE, sw=1.6)
    # middle: Python
    c.text(
        520,
        40,
        "Python's own path",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
        fill=P.PYTHON,
    )
    py = [
        "meshioplusplus.write()",
        "_helpers.write → format shim",
        "_core.<fmt>_write (pybind11)",
    ]
    for k, e in enumerate(py):
        c.box(
            520, 52 + k * 58, 200, 36, e, color=P.PYTHON, size=P.SIZE_SMALL, mono=True
        )
        if k < 2:
            c.arrow(620, 88 + k * 58, 620, 108 + k * 58, stroke=P.PYTHON, sw=1.4)
    c.arrow(620, 204, 560, 330, stroke=P.PYTHON, sw=1.6)
    # right: WASM
    c.text(
        760, 40, "WASM", size=P.SIZE_SMALL, anchor="start", weight="700", fill=P.WASM
    )
    c.box(
        760,
        52,
        176,
        36,
        "writeMesh / convert",
        color=P.WASM,
        size=P.SIZE_SMALL,
        mono=True,
    )
    c.arrow(848, 88, 848, 108, stroke=P.WASM, sw=1.4)
    c.box(
        760,
        110,
        176,
        36,
        "registry_writers() lambda",
        color=P.WASM,
        size=P.SIZE_SMALL,
        mono=True,
    )
    c.arrow(848, 146, 700, 330, stroke=P.WASM, sw=1.6)
    # the writer
    c.box(
        300,
        332,
        460,
        54,
        ["the per-format writer  (write_vtu, write_med, write_gid, …)"],
        color=P.FORMATS,
        sub="renders kProvenanceTag + the live record into the format's own comment slot (SlotTier: Block / SingleLine / Bounded / None)",
        size=P.SIZE_LABEL,
        sub_size=9.5,
    )
    # the scope
    c.box(
        24,
        300,
        236,
        120,
        ["thread-local ProvenanceScope"],
        color=P.CABI,
        sub=[
            "Mode: Off / BestEffort (default) / Required",
            "source · target · operation chain",
            "conversion-assumption notes",
            "provenance_begin_write() resets the",
            "scope-less record at each public entry",
        ],
        size=P.SIZE_SMALL,
        sub_size=9.5,
    )
    c.arrow(260, 360, 298, 360, stroke=P.CABI, sw=1.8)
    c.label(279, 348, "reads", size=9)
    # the bridge
    c.box(
        520,
        240,
        200,
        44,
        ["_provenance.scope (Python)"],
        color=P.PYTHON,
        size=P.SIZE_SMALL,
        sub="note / set_source / add_operation",
        sub_size=9,
    )
    c.arrow(520, 262, 262, 340, stroke=P.PYTHON, sw=1.4, dash="4 3", both=True)
    c.label(
        290,
        296,
        "mirrored: provenance_scope_push / _pop",
        size=9.5,
        fill=P.PYTHON,
        anchor="start",
    )
    lines = [
        "Only four write paths reach registry_write_ex, so a record carried by WriteOptions would be invisible from Python and WASM, the surfaces most callers use.",
        "A thread-local, RAII-scoped record (the set_buffer_allocator shape, per thread because writes are not serialised) is read by every writer exactly where",
        "it already read the credit line, so no signature changed. The public entry points call provenance_begin_write(), which keeps a note raised by an unrelated",
        "earlier operation out of THIS file's header; spanning several writes is what an explicit scope is for. Nothing re-emits a block read from an input file:",
        "writers render from the live record only, so provenance is replaced, never appended.",
    ]
    for k, s in enumerate(lines):
        c.label(480, 432 + 16 * k, s)
    c.legend(
        24,
        530,
        [
            (P.CORE, "C++ core surfaces"),
            (P.PYTHON, "Python"),
            (P.WASM, "WASM"),
            (P.CABI, "the record"),
            (P.FORMATS, "format writers"),
        ],
    )
    return c.render()


def interop_layers():
    c = Canvas(
        900,
        530,
        "Interoperability: a pure payload layer under thin, lazily-imported wrappers; a host-to-device move is always one bus transfer",
    )
    c.box(
        350,
        30,
        200,
        44,
        ["Mesh"],
        color=P.INK,
        sub="capsule-owned numpy from the C++ core",
        size=P.SIZE_LABEL,
        sub_size=9.5,
        fill_opacity=0.06,
    )
    c.arrow(450, 74, 450, 100, stroke=P.INK_2, sw=1.8)
    c.rrect(
        24, 102, 852, 92, fill=P.CORE, fill_opacity=0.06, stroke=P.CORE, rx=8, sw=1.4
    )
    c.text(
        40,
        124,
        "pure payload layer — imports no third-party library, mutates nothing; tested in the default CI matrix with none of the targets installed",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
        fill=P.CORE,
    )
    payloads = [
        ("_to_vtk_payload", "VTK offsets/types; region masks + sidecar"),
        ("_to_triangles_payload", "extract_surface → linearize → simplexify"),
        ("_to_table_payload", "columns; multi-component as fixed_size_list"),
        ("_to_device_payload", "block-major cells, regions as index arrays"),
    ]
    px = spread(40, 860, 4, 196)
    for x, (name, what) in zip(px, payloads):
        c.box(
            x,
            136,
            196,
            48,
            name,
            color=P.CORE,
            sub=what,
            size=P.SIZE_SMALL,
            sub_size=9,
            mono=False,
        )
    targets = [
        (
            "PyVista",
            "to_pyvista · from_pyvista",
            "shares points and data",
            P.FORMATS,
            0,
        ),
        ("trimesh", "to_trimesh · from_trimesh", "triangles only", P.FORMATS, 1),
        ("Arrow / Parquet", "to_arrow · from_arrow", "zero-copy scalars", P.DATA, 2),
        (
            "pandas · polars",
            "to_pandas · to_polars",
            "v_0/v_1 suffixes (pandas)",
            P.DATA,
            2,
        ),
        ("Blender", "to_blender · from_blender", "n-gons kept; float32", P.PYTHON, 0),
        (
            "CuPy · torch · JAX",
            "to_cupy · to_torch · to_jax",
            "one bus transfer per array",
            P.CABI,
            3,
        ),
    ]
    tx = spread(24, 876, len(targets), 136)
    for x, (name, fns, note, colour, src) in zip(tx, targets):
        c.box(
            x,
            240,
            136,
            58,
            name,
            color=colour,
            sub=[fns, note],
            size=P.SIZE_SMALL,
            sub_size=8.5,
            mono=False,
        )
        c.arrow(px[src] + 98, 184, x + 68, 238, stroke=colour, sw=1.2)
    for k, (name, note) in enumerate(
        [
            ("Open3D", "copying API, ~400 MB wheel"),
            ("DOLFINx", "one cell type, MPI, basix order"),
        ]
    ):
        c.box(
            24 + k * 176,
            322,
            166,
            40,
            name,
            color=P.MUTED,
            sub=note,
            size=P.SIZE_SMALL,
            sub_size=8.5,
            dash="4 3",
        )
    c.label(
        384,
        346,
        "Phase 2 stubs: they raise NotImplementedError by name; the constraints are recorded on the page",
        anchor="start",
    )
    c.text(
        24,
        392,
        "the zero-copy contract:",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
    )
    c.label(
        24,
        408,
        "a buffer is shared when the target takes it as-is; every copy is a note (warned), or an InteropError under zero_copy_only=True.",
        anchor="start",
    )
    c.label(
        24,
        424,
        "Derived arrays (VTK offsets/types) are not losses and never warn. Every wrapper keeps the shared arrays alive (_meshioplusplus_refs).",
        anchor="start",
    )
    c.label(
        24,
        452,
        "GPU: host → device is always a transfer across the bus; the feature removes the file round-trip and every extra copy on either side of that one transfer,",
        anchor="start",
    )
    c.label(
        24,
        468,
        "so to_cupy takes no zero_copy_only (a lie on a 100 %-copy path) and from_cupy performs one deliberate, noted device → host copy per array.",
        anchor="start",
    )
    c.label(
        24,
        484,
        "There is no [gpu], [torch] or [jax] extra: those wheels are CUDA-version-specific, so the gate names the exact pip install instead.",
        anchor="start",
    )
    c.legend(
        24,
        518,
        [
            (P.FORMATS, "mesh libraries"),
            (P.DATA, "tabular"),
            (P.PYTHON, "Blender"),
            (P.CABI, "device / tensor"),
            (P.MUTED, "not yet"),
        ],
    )
    return c.render()


def data_locations():
    c = Canvas(
        900,
        510,
        "The three data locations, and the five data operations that act on them without ever touching geometry",
    )
    # geometry
    c.box(
        24,
        50,
        180,
        110,
        ["geometry"],
        color=P.INK,
        sub=["points · cells · regions", "never modified by a data op"],
        size=P.SIZE_LABEL,
        sub_size=9.5,
        fill_opacity=0.05,
    )
    # shelves
    shelves = [
        ("point_data", "one array per point: (N,) or (N, k)", P.DATA),
        ("cell_data", "a LIST: one array per cell block, aligned to cells", P.DATA),
        ("field_data", "whole-mesh values (scalars, small arrays)", P.DATA),
    ]
    for k, (name, what, colour) in enumerate(shelves):
        y = 50 + k * 60
        c.box(
            250,
            y,
            300,
            46,
            name,
            color=colour,
            sub=what,
            size=P.SIZE_LABEL,
            sub_size=9.5,
            mono=False,
        )
    c.line(204, 105, 250, 73, stroke=P.MUTED, sw=1, dash="3 3")
    c.line(204, 105, 250, 133, stroke=P.MUTED, sw=1, dash="3 3")
    c.label(227, 120, "aligned", size=9, rotate=-35)
    # ops
    ops = [
        (
            "data_manage",
            "keep → drop → rename",
            "values pass through byte-identical",
            P.CORE,
            [0, 1, 2],
        ),
        (
            ["point_data_to_cell_data", "cell_data_to_point_data"],
            "averaging; always Float64",
            "serial scatter, so thread-count-independent",
            P.CABI,
            [0, 1],
        ),
        (
            "data_calc",
            "expression → new array",
            "own parser: + − * / abs sqrt min max norm",
            P.FORMATS,
            [0, 1],
        ),
        (
            "data_condition",
            "clamp / normalize / standardize",
            "statistics joint across blocks",
            P.PYTHON,
            [0, 1],
        ),
        (
            "data_info",
            "read-only report",
            "never throws on NaN/inf: it counts them",
            P.MUTED,
            [0, 1, 2],
        ),
    ]
    for k, (name, what, note, colour, locs) in enumerate(ops):
        y = 50 + k * 60
        c.box(
            600,
            y,
            276,
            52,
            name,
            color=colour,
            sub=[what, note],
            size=P.SIZE_SMALL,
            sub_size=8.5,
            mono=False,
        )
        for loc in locs:
            c.line(
                550, 50 + loc * 60 + 23, 600, y + 26, stroke=colour, sw=1, opacity=0.6
            )
    c.text(
        24, 368, "cross-cutting rules", size=P.SIZE_SMALL, anchor="start", weight="700"
    )
    rules = [
        "non-finite values are always excluded from every reduction; NanPolicy (ignore / replace / fail) only decides what reaches the output",
        "an op writing cell_data emits exactly NumCellBlocks() arrays; an op reading one validates the block count first",
        "components = product of the trailing dims; 1-D scalars stay 1-D; point/cell sets never enter the core",
    ]
    for k, s in enumerate(rules):
        c.label(24, 386 + 16 * k, "• " + s, anchor="start")
    c.text(
        24,
        446,
        "not data operations, though the CLI groups them under `data`:",
        size=P.SIZE_SMALL,
        anchor="start",
        weight="700",
    )
    c.label(
        24,
        462,
        "gradient · hessian · estimate_error · data_integrate read geometry (face areas, cell measures, adjacency): mesh operations that happen to produce arrays",
        anchor="start",
    )
    c.label(
        24,
        478,
        "names may contain ':' (gmsh:physical, refine:level); the CLI splits `data rename OLD:NEW` on the LAST colon and `data calc \"NAME = EXPR\"` on the FIRST '='",
        anchor="start",
    )
    c.legend(
        24,
        496,
        [
            (P.INK, "geometry"),
            (P.DATA, "data locations"),
            (P.CORE, "manage"),
            (P.CABI, "average"),
            (P.FORMATS, "calc"),
            (P.PYTHON, "condition"),
            (P.MUTED, "info"),
        ],
    )
    return c.render()


def ml_fanout():
    c = Canvas(
        900,
        420,
        "From meshes to training tensors: graphs, feature matrices, datasets and framework tensors",
    )
    sources = [
        ("Mesh", "one solution"),
        ("TimeSeries", "a transient run; nothing cached"),
        ("DatasetManifest", "many runs, splits, tags; JSON"),
    ]
    for k, (name, what) in enumerate(sources):
        c.box(
            24,
            50 + k * 70,
            190,
            54,
            name,
            color=P.INK,
            sub=what,
            size=P.SIZE_LABEL,
            sub_size=9,
            fill_opacity=0.05,
            mono=False,
        )
    mids = [
        (
            "edge_index(kind, undirected)",
            "(2, E) int64 — node graph from _REFINE_EDGES, or the cell dual",
            P.CORE,
        ),
        (
            "feature_matrix(location, fields, regions)",
            "(N, F) float64 + a versioned column schema",
            P.DATA,
        ),
        (
            "write_dataset(source, path, format)",
            "parquet (hive) · zarr · hdf5; strict schema; manifest last",
            P.FORMATS,
        ),
        (
            "to_torch · to_jax · to_dlpack",
            "framework tensors; host adoption zero-copy via DLPack",
            P.CABI,
        ),
    ]
    for k, (name, what, colour) in enumerate(mids):
        y = 40 + k * 66
        c.box(
            290,
            y,
            320,
            54,
            name,
            color=colour,
            sub=what,
            size=P.SIZE_SMALL,
            sub_size=9,
            mono=False,
        )
        for s in range(3):
            c.line(214, 77 + s * 70, 290, y + 27, stroke=colour, sw=0.9, opacity=0.5)
    outs = [
        (
            "PyG Data / DataLoader",
            "graph_sample: pos, x, y, edge_index, edge_attr",
            P.CORE,
        ),
        ("PhysicsNeMo Reader", "make_reader / make_dataset; stats JSON", P.WASM),
        (
            "to_physicsnemo / from_physicsnemo",
            "physicsnemo.mesh.Mesh (one simplex kind)",
            P.WASM,
        ),
    ]
    for k, (name, what, colour) in enumerate(outs):
        y = 60 + k * 80
        c.box(
            660,
            y,
            216,
            60,
            name,
            color=colour,
            sub=what,
            size=P.SIZE_SMALL,
            sub_size=8.5,
            mono=False,
        )
    c.arrow(610, 67, 660, 90, stroke=P.CORE, sw=1.3)
    c.arrow(610, 133, 660, 100, stroke=P.DATA, sw=1.3)
    c.arrow(610, 199, 660, 170, stroke=P.FORMATS, sw=1.3)
    c.arrow(610, 265, 660, 250, stroke=P.CABI, sw=1.3)
    c.label(
        450,
        322,
        "the streaming invariant holds everywhere: write_dataset and the readers use the same per-entry read loop as read_sequence,",
    )
    c.label(
        450,
        338,
        "so one mesh is alive at a time (two for a t → t+n pair); the feature-matrix column order is the versioned contract (FEATURE_SCHEMA_VERSION):",
    )
    c.label(
        450,
        354,
        "coords → data arrays (sorted, or the explicit order) → region one-hots; a schema mismatch between entries is an error, never a silent union",
    )
    c.legend(
        24,
        400,
        [
            (P.CORE, "graph"),
            (P.DATA, "features"),
            (P.FORMATS, "datasets on disk"),
            (P.CABI, "tensors"),
            (P.WASM, "PhysicsNeMo"),
        ],
    )
    return c.render()


def abi_tiers():
    c = Canvas(
        760,
        590,
        "Which C++ header changes bump MESHIOPLUSPLUS_ABI_VERSION: the Tier A / B / C decision",
    )

    def diamond(cx, cy, w, h, lines, colour):
        c.polygon(
            [(cx, cy - h / 2), (cx + w / 2, cy), (cx, cy + h / 2), (cx - w / 2, cy)],
            fill=colour,
            fill_opacity=0.10,
            stroke=colour,
            sw=1.5,
        )
        for k, s in enumerate(lines):
            c.text(
                cx,
                cy - (len(lines) - 1) * 7 + k * 14 + 4,
                s,
                size=P.SIZE_SMALL,
                weight="600" if k == 0 else None,
                fill=P.INK if k == 0 else P.INK_2,
            )

    c.box(
        230,
        30,
        300,
        44,
        ["a header under src/cpp/include/meshioplusplus/ changed"],
        color=P.INK,
        sub="all of them are installed and compiled by consumers; src/cpp/src/**.cpp never matters",
        size=P.SIZE_SMALL,
        sub_size=9,
        fill_opacity=0.05,
    )
    c.arrow(380, 74, 380, 100, stroke=P.INK_2)
    diamond(
        380,
        160,
        360,
        116,
        [
            "does it change a type's LAYOUT?",
            "data member added / removed / retyped,",
            "base or virtual, alignment,",
            "enum underlying type, template params",
        ],
        P.ERROR,
    )
    c.text(575, 150, "yes", size=P.SIZE_SMALL, weight="700", fill=P.ERROR)
    c.arrow(560, 160, 600, 160, stroke=P.ERROR)
    c.box(
        602,
        128,
        140,
        64,
        ["Tier A", "bump"],
        color=P.ERROR,
        sub="test_abi_layout.cpp catches it",
        sub_size=8.5,
    )
    c.text(392, 236, "no", size=P.SIZE_SMALL, weight="700", fill=P.INK_2)
    c.arrow(380, 218, 380, 254, stroke=P.INK_2)
    diamond(
        380,
        316,
        360,
        116,
        [
            "does it edit the BODY of something",
            "that already existed?",
            "an inline function, a template,",
            "a default argument",
        ],
        P.PYTHON,
    )
    c.text(575, 306, "yes", size=P.SIZE_SMALL, weight="700", fill=P.PYTHON)
    c.arrow(560, 316, 600, 316, stroke=P.PYTHON)
    c.box(
        602,
        284,
        140,
        64,
        ["Tier B", "bump"],
        color=P.PYTHON,
        sub="ODR: only a recorded review",
        sub_size=8.5,
    )
    c.text(392, 392, "no", size=P.SIZE_SMALL, weight="700", fill=P.INK_2)
    c.arrow(380, 374, 380, 410, stroke=P.INK_2)
    c.box(
        230,
        412,
        300,
        64,
        ["Tier C — additive, do NOT bump"],
        color=P.WASM,
        sub=[
            "a new inline function, constexpr, type, declaration or header;",
            "an appended enumerator (CellType has a fixed underlying type)",
        ],
        sub_size=8.5,
    )
    c.arrow(530, 444, 600, 444, stroke=P.WASM)
    c.box(
        602,
        412,
        140,
        64,
        ["record it in", "abi_reviews.md"],
        color=P.WASM,
        sub="release + every changed header",
        sub_size=8.5,
    )
    lines = [
        "tools/check-abi-version.sh (the lint job) fails a build whose installed headers changed while the ABI version did not",
        "and no review row names that release and those headers. The ABI version is deliberately NOT the release version:",
        "it moves only when already-compiled consumers stop being compatible (11 at v10.21.0), and it is the C++ variants' SOVERSION.",
        "The C API is a separate promise: SOVERSION 0, SameMajorVersion, append-only option structs with reserved tails,",
        "so a C consumer pins find_package(meshioplusplus 10 … COMPONENTS C) and a C++ consumer pins the exact release or the ABI number.",
    ]
    for k, s_ in enumerate(lines):
        c.label(380, 504 + 16 * k, s_)
    return c.render()


def roadmap_map():
    c = Canvas(
        960,
        560,
        "The roadmap at a glance: open items by theme, shaded by effort, with their dependencies",
    )
    columns = [
        (
            "§1 dashboard + training",
            [
                ("overview · drill-down · health", "M"),
                ("companion service — design first", "L"),
                ("launch · monitor · history · logs", "M"),
                ("prediction preview in-viewer", "M"),
            ],
            [(1, 2)],
        ),
        (
            "§2 scale",
            [
                ("10M+ cell benchmark tier", "S"),
                ("streaming writes of one mesh", "L"),
                ("out-of-core operations", "XL"),
            ],
            [(0, 1), (1, 2)],
        ),
        (
            "§3 ecosystem",
            [
                ("Rust bindings over the C API", "M"),
                ("conda-forge · CRAN · Julia General", "M"),
                ("Blender / ParaView listings", "S"),
            ],
            [],
        ),
        (
            "§4 quality",
            [
                ("fuzzing the 43 readers", "M"),
                ("format conformance matrix", "M"),
                ("property-based tests", "M"),
            ],
            [],
        ),
        (
            "§5 NURBS — spike first",
            [
                ("spike: can the model stretch?", "M"),
                ("read-only CAD ingestion", "L"),
                ("a real IGA data model", "XL"),
            ],
            [(0, 1), (1, 2)],
        ),
        (
            "§6 generation",
            [
                ("box · sphere · cylinder · disk", "S"),
                ("extrude", "M"),
                ("revolve", "M"),
                ("constrained 2-D meshing", "L"),
            ],
            [(0, 1), (1, 2)],
        ),
        ("§7 chat", [("chat verb over TOOL_REGISTRY", "M")], []),
        (
            "§8 binding parity",
            [
                ("decimate_volume bindings", "S"),
                ("flat-ABI gaps (frozen, sets, …)", "M"),
            ],
            [],
        ),
        (
            "§9 formats",
            [
                ("exodus writer: sets, multi-step", "M"),
                ("regions Phase 2 formats", "M"),
                ("read_med_metadata", "S"),
                ("gmsh $Periodic / 4.0, pyramid14", "M"),
            ],
            [],
        ),
        (
            "§10 follow-ups elsewhere",
            [("pipeline · gpu · interop · julia · mcp", "S")],
            [],
        ),
    ]
    per_row = 5
    col_w = 176
    col_gap = (960 - 48 - per_row * col_w) / (per_row - 1)
    box_h = 40
    for k, (title, items, deps) in enumerate(columns):
        row, col = divmod(k, per_row)
        x = 24 + col * (col_w + col_gap)
        y0 = 40 + row * 250
        c.text(x, y0, title, size=P.SIZE_SMALL, anchor="start", weight="700")
        centres = []
        for i, (label, effort) in enumerate(items):
            y = y0 + 10 + i * (box_h + 10)
            colour = P.EFFORT[effort]
            dashed = "design first" in label or "spike" in label
            c.rrect(
                x,
                y,
                col_w,
                box_h,
                fill=colour,
                fill_opacity=0.9 if effort in ("L", "XL") else 0.55,
                stroke=colour,
                rx=5,
                sw=1.2,
                dash="5 3" if dashed else None,
            )
            c.text(
                x + 6,
                y + 17,
                label,
                size=9.5,
                anchor="start",
                fill=P.WHITE if effort in ("L", "XL") else P.INK,
                weight="600",
            )
            c.text(
                x + 6,
                y + 31,
                effort,
                size=9,
                anchor="start",
                fill=P.WHITE if effort in ("L", "XL") else P.INK_2,
                mono=True,
            )
            centres.append((x + col_w - 12, y + box_h))
        for a, b in deps:
            ax, ay = centres[a]
            c.arrow(ax, ay, ax, ay + 10, stroke=P.INK, sw=1.4)
    c.legend(
        24,
        526,
        [
            (P.EFFORT["S"], "S: days"),
            (P.EFFORT["M"], "M: a couple of weeks"),
            (P.EFFORT["L"], "L: a month or more"),
            (P.EFFORT["XL"], "XL: a project in its own right"),
        ],
    )
    c.label(
        24,
        546,
        "dashed border: a design pass or research spike must precede the item; an arrow means the lower box depends on the one above it",
        anchor="start",
    )
    return c.render()


FIGURES = {
    "pipeline_flow": pipeline_flow,
    "sequences_shapes": sequences_shapes,
    "provenance_write_paths": provenance_write_paths,
    "interop_layers": interop_layers,
    "data_locations": data_locations,
    "ml_fanout": ml_fanout,
    "abi_tiers": abi_tiers,
    "roadmap_map": roadmap_map,
}
