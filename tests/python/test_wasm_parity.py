"""Cross-checks between the WASM binding's four textual surfaces.

@meshioplusplus/wasm has no build step in the default CI matrix (it needs
Emscripten), so nothing catches its own internal drift: the raw embind
registration (``bindings/wasm/js_bindings.cpp``), the ergonomic wrapper
(``src/wasm/src/index.mjs``), its ambient types (``src/wasm/index.d.ts``) and
the C API (``bindings/c/include/meshioplusplus/meshioplusplus.h``) are four
independent, hand-maintained texts that are supposed to describe the same
surface. This is pure text parsing (regexes over the four files) so it runs
in the normal Python job like any other test.

See doc/roadmap.md (closed WASM-parity item) and AGENTS.md's ABI/bindings
section for the broader convention this enforces.
"""

import pathlib
import re

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
JS_BINDINGS = REPO / "bindings" / "wasm" / "js_bindings.cpp"
INDEX_MJS = REPO / "src" / "wasm" / "src" / "index.mjs"
INDEX_DTS = REPO / "src" / "wasm" / "index.d.ts"
C_HEADER = REPO / "bindings" / "c" / "include" / "meshioplusplus" / "meshioplusplus.h"

# The two stateful handles: the transient-XDMF writer (9 free functions
# folded into `createXdmfTimeSeriesWriter`'s returned object) and the
# sequence reader (8 free functions folded into `openSequence`'s returned
# object) -- both embind-registered but deliberately with no top-level
# wrapper key of their own name.
_RAW_ONLY = {
    "xdmfSeriesCreate",
    "xdmfSeriesWritePointsCells",
    "xdmfSeriesWriteData",
    "xdmfSeriesWriteDataArrays",
    "xdmfSeriesFlush",
    "xdmfSeriesFinalize",
    "xdmfSeriesNumSteps",
    "xdmfSeriesFinalized",
    "xdmfSeriesFree",
    "sequenceOpen",
    "sequenceCount",
    "sequencePath",
    "sequenceStep",
    "sequenceTime",
    "sequenceTimeSource",
    "sequenceRead",
    "sequenceFree",
}

# Wrapper keys with no `mio_*` counterpart, and why.
_JS_ONLY = {
    "FS": "Emscripten virtual-filesystem passthrough, not a mesh operation",
    "withProvenance": "JS-only composition of provenanceBegin/provenanceEnd",
    "availableFormats": "bulk introspection; C exposes mio_format_readable/writable per-format",
    "hasCgnslib": "build-capability introspection helper",
    "numNodesPerCell": "bulk map; C exposes mio_cell_type_num_nodes per cell type",
    "parallelBackend": "build-capability introspection helper",
    "topologicalDimension": "bulk map; C exposes mio_cell_type_dimension per cell type",
    "convertSurface": "WASM/browser rendering convenience, no C entry point",
    "convertSurfaceOps": "WASM/browser rendering convenience, no C entry point",
    "createXdmfTimeSeriesWriter": "composes mio_xdmf_series_* (see _RAW_ONLY)",
    "readProvenance": "shaped differently on C: mio_read_metadata_provenance_*",
    "sequenceEntries": "bulk listing convenience; C's sequence API is the stateful open/count/path/... handle",
    "openSequence": "composes mio_sequence_* (see _RAW_ONLY)",
}

# camelCase wrapper key -> mio_<name> (without the mio_ prefix / _ex suffix,
# which the check adds itself) when it is not the mechanical snake_case of
# the wrapper key.
_JS_TO_C = {
    "readMesh": "read",
    "readMeshSelective": "read",
    "readMetadata": "read_metadata_create",
    "dataInfo": "data_info_create",
    "dataIntegrate": "data_integrate_create",
    "writeMesh": "write",
    "provenanceBegin": "provenance_scope_begin",
    "provenanceEnd": "provenance_scope_end",
    "runPipeline": "pipeline_run_json",
}

# `mio_*` operation entry points (first parameter `const mio_mesh*`) with no
# WASM binding, and why.
_C_WITHOUT_JS = {
    # Aggregate scalar counts with no dedicated JS need: `attachQuality`
    # already exposes the per-cell metrics as cell_data, and
    # `surfaceWatertightCheck` covers the surface-quality summary case.
    "quality_counts",
}


def _read(path):
    return path.read_text(encoding="utf-8")


def _embind_names():
    return re.findall(r'emscripten::function\("(\w+)"', _read(JS_BINDINGS))


def _wrapper_keys():
    text = _read(INDEX_MJS)
    start = text.index("export async function loadMeshioPlusPlus")
    body = text[start:]
    ret_start = body.index("return {")
    ret_end = body.index("\n    };", ret_start)
    block = body[ret_start:ret_end]
    return re.findall(r"^ {8}(\w+):", block, re.MULTILINE)


def _jsdoc_returns_names():
    text = _read(INDEX_MJS)
    start = text.index("@returns {Promise<{")
    end = text.index(" * }>}", start)
    block = text[start:end]
    return re.findall(r"^\s*\*\s+(\w+)(?:<\w+>)?:", block, re.MULTILINE)


def _dts_member_names():
    text = _read(INDEX_DTS)
    start = text.index("export interface MeshioPlusPlusModule {")
    depth = 0
    end = start
    for i, ch in enumerate(text[start:], start):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                end = i
                break
    block = text[start:end]
    return re.findall(r"^ {2}(\w+)(?:<\w+>)?[(:]", block, re.MULTILINE)


def _c_declarations():
    text = _read(C_HEADER)
    pat = re.compile(r"MIO_API\s+[\w*\s]+?\bmio_([a-zA-Z0-9_]+)\s*\(\s*([^,)]*)")
    seen = {}
    for m in pat.finditer(text):
        name, first_param = m.group(1), m.group(2).strip()
        seen.setdefault(name, first_param)
    return seen


def _c_operation_names(declarations):
    """Names of `mio_*` entry points shaped like a mesh operation: first
    parameter is a `mio_mesh*`, excluding mesh accessors/mutators
    (`mio_mesh_*`) and result-object constructors (`*_create`)."""
    return {
        name
        for name, first_param in declarations.items()
        if re.search(r"\bmio_mesh\s*\*", first_param)
        and not name.startswith("mesh_")
        and not name.endswith("_create")
    }


def _camel_to_snake(name):
    return re.sub(r"(?<!^)(?=[A-Z])", "_", name).lower()


def _snake_to_camel(name):
    head, *tail = name.split("_")
    return head + "".join(p.title() for p in tail)


def test_embind_functions_are_reachable_from_js():
    embind = set(_embind_names())
    wrapper = set(_wrapper_keys())
    unreachable = embind - wrapper - _RAW_ONLY
    assert not unreachable, (
        f"embind function(s) with no wrapper key and not in _RAW_ONLY: "
        f"{sorted(unreachable)}"
    )


def test_raw_only_exemptions_are_not_stale():
    embind = set(_embind_names())
    stale = _RAW_ONLY - embind
    assert not stale, f"_RAW_ONLY names no longer embind-registered: {sorted(stale)}"


def test_wrapper_dts_and_jsdoc_agree():
    wrapper = set(_wrapper_keys())
    dts = set(_dts_member_names())
    jsdoc = set(_jsdoc_returns_names())
    assert wrapper, "failed to extract any wrapper keys from index.mjs"
    assert dts, "failed to extract any interface members from index.d.ts"
    assert jsdoc, "failed to extract any @returns typedef names from index.mjs"
    assert wrapper == dts, f"wrapper vs index.d.ts mismatch: {wrapper ^ dts}"
    assert wrapper == jsdoc, f"wrapper vs JSDoc @returns mismatch: {wrapper ^ jsdoc}"


def test_wrapper_keys_map_to_c_entry_points():
    declarations = _c_declarations()
    assert declarations, "failed to extract any mio_* declarations"
    for key in _wrapper_keys():
        if key in _JS_ONLY:
            continue
        base = _JS_TO_C.get(key, _camel_to_snake(key))
        assert base in declarations or f"{base}_ex" in declarations, (
            f"wrapper key '{key}' has no mio_{base}/mio_{base}_ex C entry point; "
            "add one, fix _JS_TO_C, or add a reasoned _JS_ONLY exemption"
        )


def test_js_only_exemptions_are_not_stale():
    declarations = _c_declarations()
    wrapper = set(_wrapper_keys())
    for key in _JS_ONLY:
        assert key in wrapper, f"_JS_ONLY names a key no longer in the wrapper: {key}"
        base = _JS_TO_C.get(key, _camel_to_snake(key))
        assert (
            base not in declarations and f"{base}_ex" not in declarations
        ), f"_JS_ONLY['{key}'] is stale: mio_{base} now exists"


def test_c_operations_are_reachable_from_js():
    declarations = _c_declarations()
    wrapper = set(_wrapper_keys())
    c_to_js = {v: k for k, v in _JS_TO_C.items()}
    for name in sorted(_c_operation_names(declarations)):
        if name in _C_WITHOUT_JS:
            continue
        # A `_ex` (opts-struct) entry point is the same operation as its
        # base name for JS purposes: the wrapper always takes the full
        # parameter set, so there is no separate "Ex" wrapper key.
        base_name = name[: -len("_ex")] if name.endswith("_ex") else name
        camel = c_to_js.get(base_name, _snake_to_camel(base_name))
        assert camel in wrapper, (
            f"mio_{name} has no WASM wrapper entry (expected '{camel}'); "
            "add one or add a reasoned _C_WITHOUT_JS exemption"
        )


def test_c_without_js_exemptions_are_not_stale():
    declarations = _c_declarations()
    operations = _c_operation_names(declarations)
    stale = _C_WITHOUT_JS - operations
    assert not stale, f"_C_WITHOUT_JS names no longer operation-shaped: {sorted(stale)}"
