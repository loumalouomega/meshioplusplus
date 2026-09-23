# FEBio plot file (`.xplt`)

The `.xplt` file is the binary results database [FEBio](https://febio.org) writes during a run and FEBio Studio opens. It holds the mesh once and then one **state** per saved time. meshio++ reads it one state at a time, so a run converts to a `.vtu` sequence. The model's input is the separate [`.feb`](./febio.md) file.

| | |
|---|---|
| **Format name** | `xplt` |
| **Extensions** | `.xplt` (also recognised by content: its magic, in either byte order) |
| **Read / Write** | ✓ / — (read-only) |
| **Extra dependencies** | zlib for compressed files: built into the C++ core with `-DMESHIOPLUSPLUS_WITH_ZLIB=ON`; the Python fallback uses the standard library's |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("run.xplt")                  # the first state
last = meshioplusplus.read("run.xplt", time_step=-1)    # the last one
meshioplusplus.xplt.time_values("run.xplt")             # every state's time
for time, mesh in meshioplusplus.read_sequence("run.xplt"):
    ...
```

```bash
meshioplusplus convert run.xplt 'out_{step}.vtu'   # one .vtu per state
meshioplusplus info run.xplt                       # blocks, regions, arrays, times
```

`read` also takes `points_only`/`arrays` (read only the named arrays) and `lenient` (read FEBio's `tet5`/`tet15` domains as `tetra`/`tetra10`). Both engines read every state alike.

## File structure

After a 4-byte magic (`0x00464542`), the file is a tree of chunks, each a `u32` id, a `u32` payload size and the payload. The layout follows FEBio's writer (`FEBioPlot/FEBioPlotFile.cpp`, `PltArchive.cpp`) and FEBio Studio's reader (`XPLTLib/xpltReader3.cpp`):

- a root with the header (plot version, compression flag, software) and the **dictionary**: every plot variable's name, type (`float`, `vec3f`, `mat3fs`, `mat3fd`, `mat3f`, `tens4fs`, arrays) and storage format;
- the **mesh**: nodes (from version 0x0033 with their ids), domains (element type, part, elements), surfaces, node sets, element sets, facet sets, parts;
- one **state** per saved time: the time, a status, and each variable's values per region.

When compression is on, every chunk after the first mesh (every state) is its own zlib stream. Plot versions `0x0030` and later are read, which covers FEBio 3 and 4 (FEBio 4.12 writes `0x0035`), in either byte order. FEBio 2's plot files are refused, naming the version.

## Mapping

| `.xplt` | meshio++ |
|---|---|
| a domain | a cell block, and a `cell` region tagged with its part id. The region is named after the domain; FEBio 4 names only solid domains, so a nameless one takes the name of the element set holding exactly its elements (FEBio writes one per `<Elements>` block), else its part's. |
| node set / element set | `point` / `cell` region |
| surface, facet set | a `side` region when every facet is a face of a solid, otherwise its own cell block and a `cell` region, as for [`.feb`](./febio.md#mapping) |
| the state's time, index, status | `field_data["meshio:time"]`, `"xplt:step"`, `"xplt:status"` |
| nodal variables | point data |
| domain variables stored per element (`FMT_ITEM`) | cell data, NaN on domains that do not carry the variable |
| domain variables stored once per domain (`FMT_REGION`) | cell data, the value repeated over the domain's cells |
| domain variables stored per element node (`FMT_NODE`, `FMT_MULT`) | point data, averaged over the elements that share each node |
| global variables | field data |
| surface, edge and material-point variables | not read; one warning names them |

Values are FEBio's `float32`, widened to `float64`. Components keep FEBio's order. A symmetric tensor (`mat3fs`) has 6 components, `xx, yy, zz, xy, yz, xz`, which is also VTK's order. A diagonal tensor has 3 components, and a full one (`mat3f`) 9 in row-major order.

**Compared with FEBio Studio's VTK export.** FEBio Studio exports node data and per-element data the same way. For per-element-node data it keeps the value of the last element written at a shared node, where meshio++ averages; the two agree for continuous fields. FEBio Studio drops `mat3f` and `tens4fs` data, while meshio++ keeps them. It writes `hex27` as a 20-node hexahedron, while meshio++ keeps all 27 nodes through the `"febio"` [node-ordering](../node_ordering.md) table.

## Notes

- **A state cut short** (FEBio killed mid-write) is dropped with a warning, and the earlier states are still read.
- **A remeshed run**, one with a second mesh between states, is a `ReadError`.
- **Validation**, made outside the repository with FEBio 4.12 built from source: the displacement read from FEBio's plot files equals FEBio's own `node_data` log to `float32` precision at every state, for three models. febio-python's reader (MIT) gives identical displacements, strains, stresses and fibre vectors for its `0x0031` samples. The test fixtures are plot files FEBio 4.12 wrote for the models `tools/gen_febio_fixtures.py` generates, one of them compressed.
