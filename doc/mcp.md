# MCP server

Expose every meshio++ operation to AI agents over the
[Model Context Protocol](https://modelcontextprotocol.io/): reading and
writing 40+ mesh formats, conversion, and the full mesh- and data-operation
suite become **tools** any MCP client (Claude Code, Claude Desktop, the MCP
inspector, …) can call.

```bash
pip install "meshioplusplus[mcp]"          # the mcp SDK needs Python >= 3.10
claude mcp add meshioplusplus -- meshioplusplus-mcp
```

Then ask the agent things like *"convert `bracket.msh` to VTU, report its
quality, and slice it at z = 0.02"* — it drives `convert`, `quality` and
`slice` itself.

Every tool is **stateless and file-path based**: input path(s) in, output
path(s) out, a strict-JSON report back. That mirrors the CLI, keeps arbitrarily
large meshes out of the protocol, and lets the agent work in its own filesystem
workspace. Nothing here is part of the C++ core, which stays dependency-free.

## Installation

| Extra | Brings | For |
|---|---|---|
| `meshioplusplus[mcp]` | the official [MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk) (MIT) | the `meshioplusplus-mcp` server |

This extra is **not** in `meshioplusplus[all]` — `[all]` means "the optional
dependencies the *formats* need", the same reasoning as `[interop]`/`[viewer]`.
Note the `mcp` SDK itself requires **Python ≥ 3.10** while meshio++ supports
3.8: the package (and the server's pure tool layer) works everywhere, only
*running* the server needs the newer Python.

Without the extra, everything degrades by name rather than by traceback:

```
$ meshioplusplus-mcp
meshio++: mcp server: mcp is not installed; install it with
`pip install meshioplusplus[mcp]` (requires Python >= 3.10)
```

`meshioplusplus.mcp.has_mcp()` answers the same question without raising.

## Running the server

The server speaks stdio (the standard local-agent transport):

```bash
meshioplusplus-mcp                          # unrestricted paths
meshioplusplus-mcp --root /path/to/work     # sandboxed (recommended)
python -m meshioplusplus.mcp                # equivalent
```

Claude Desktop (`claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "meshioplusplus": {
      "command": "meshioplusplus-mcp",
      "args": ["--root", "/path/to/workspace"]
    }
  }
}
```

Interactive browsing: `npx @modelcontextprotocol/inspector meshioplusplus-mcp`.

## Path sandbox

By default paths are unrestricted — the server runs locally under your own
account and MCP clients gate filesystem access themselves. Pass `--root DIR`
(or set `MESHIOPLUSPLUS_MCP_ROOT`) to confine every input **and** output path:
each path is realpath-resolved (symlinks included) and must stay inside the
root; relative paths resolve against it. Violations come back as clean
`{"error": ..., "error_type": "ValueError"}` payloads the agent can act on —
no tool ever surfaces a raw traceback.

## Tools

33 tools; the two marked *gated* need a further extra and return a named
install error without it. Transforming tools take `input_path`/`output_path`
(+ optional `input_format`/`output_format`, otherwise inferred from the
extension) and return the written path plus a mesh summary and the operation's
report.

### Inspection (read-only)

| Tool | Returns |
|---|---|
| `formats` | readable/writable format lists + extension map (also the `meshioplusplus://formats` resource) |
| `sniff` | format identified from leading bytes + extension |
| `info` | fast file summary via `read_metadata` — counts, cell blocks, data names, regions, time steps |
| `stats` | bbox, centroid, areas/volumes, per-type counts, inverted cells |
| `quality` | per-metric summaries + histograms; pass `output_path` to write the mesh with `quality:<metric>` cell data |
| `data_info` | every data array's dtype/shape/ranges/NaN counts |
| `regions` | named point/cell/side groups with kind/dim/tag and an entries preview |
| `bandwidth` | node-numbering bandwidth |
| `data_preview` | a bounded window (`offset`/`limit`) of one data array's values |
| `diff` | verdict + `equal` boolean + per-section detail for two files |

### Conversion

| Tool | Notes |
|---|---|
| `convert` | any-to-any format conversion; `points_only`/`arrays`/`time_step` narrow the read; `mode: ascii\|binary` and `compression: zlib\|lz4\|zstd\|lzma\|gzip\|none` subsume the CLI's `ascii`/`binary`/`compress`/`decompress` verbs |

### Mesh operations

`extract_surface`, `extract_skin`, `reorder` (reports bandwidth before/after),
`clean`, `crop` (bbox or half-space), `slice`, `isosurface`, `transform`,
`convert_cells`, `refine` (uniform, or a subset via `cells`/`region`/`where`
with a conforming `closure`), `decimate`, `smooth`, `merge` (N inputs), `split`
(one file per piece, `name_template`), `partition` (one file per part),
`interpolate` (source → target field transfer), `gradient` (the gradient,
divergence or curl of a `point_data` field — see
[field derivatives](/gradient); reports `num_skipped` and `num_fallback`).
Parameters mirror the Python API / CLI one-to-one; operations that produce
reports (`clean`, `decimate`, `smooth`, `gradient`) include them in the
response.

### Data operations

| Tool | Notes |
|---|---|
| `data_manage` | keep/drop/rename arrays: `keep`/`drop` are `[location, name]` pairs, `rename` is `[location, old, new]` triples |
| `data_convert` | average between locations (`direction: point_to_cell \| cell_to_point`) |
| `data_calc` | expression evaluator; accepts the CLI's `"NAME = EXPR"` spelling |
| `data_condition` | clamp / normalize / standardize |

### Gated

| Tool | Extra | Notes |
|---|---|---|
| `data_export` | `[arrow]` | data arrays → Parquet table |
| `screenshot` | `[viewer]` | off-screen PNG render, returned as MCP image content |

## Reports are strict JSON

Every response survives `json.dumps(..., allow_nan=False)`: numpy scalars and
arrays are converted, `NaN`/`±Inf` become `null` (with a
`non_finite_replaced` count so the loss is visible), and any array longer
than 1000 elements is replaced by a
`{"truncated": true, "size", "shape", "dtype", "preview"}` wrapper — reports
stay agent-sized no matter how large the mesh.

## Architecture: the pure payload layer

`src/python/meshioplusplus/mcp/` is split exactly the way
[`_interop.py`](./interop) is, and for the same reason. Everything a tool
*does* lives in the pure layer `_tools.py` — imports only meshioplusplus +
numpy + stdlib, runs on every supported Python, tested in the default CI
matrix with the `mcp` SDK absent. The FastMCP layer `_server.py` contributes
only typed signatures (the JSON schemas) and docstrings (the tool
descriptions), and is the one module importing the SDK.

`_tools.TOOL_REGISTRY` is the single source of truth: the server registers
from it, and the **parity guard** in `tests/python/test_mcp.py` asserts every
public operation in `meshioplusplus.__all__` is claimed by some tool — a new
operation fails CI until it gets a tool (or a conscious exemption). That test
is what keeps this page's tool table honest as the library grows.
