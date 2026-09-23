# FEBio input (`.feb`)

The `.feb` file is the input of [FEBio](https://febio.org), the finite-element solver for biomechanics, and what FEBio Studio saves. It is XML with a `febio_spec` root, and it holds a whole model: materials, boundary conditions, loads and steps as well as the mesh. meshio++ reads and writes **the mesh**. FEBio's results are the separate [`.xplt`](./xplt.md) plot file.

| | |
|---|---|
| **Format name** | `febio` |
| **Extensions** | `.feb` (also recognised by content: a `<febio_spec` root) |
| **Read / Write** | ✓ (spec 2.5, 3.0, 4.0) / ✓ (spec 4.0) |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.feb")                  # the mesh, sets, surfaces and MeshData
mesh = meshioplusplus.febio.read("old.feb", lenient=True)  # downgrade tet5/tet15 instead of failing
meshioplusplus.write("mesh.feb", mesh)                   # a spec-4.0 mesh
```

Both engines (the C++ core and the pure-Python reference) read the same meshes and write the same bytes.

## Spec versions

The mesh lives in `<Geometry>` in spec 2.5 and in `<Mesh>`, plus `<MeshDomains>`, in 3.0 and 4.0. The reader dispatches on the `version` attribute and follows the rules of FEBio's own parsers (`FEBioXML/FEBioGeometrySection.cpp`, `FEBioMeshSection.cpp`, `FEBioMeshSection4.cpp`). Any other version is a `ReadError` naming the three it accepts. A `<Mesh from="other.feb"/>` section is read from the other file, as FEBio does; a relative path is taken relative to the `.feb`. The spec-2.5/3.0 `<Part>`/`<Instance>` form is refused.

## Mapping

| `.feb` | meshio++ |
|---|---|
| `<Nodes>` (one or more; ids may be sparse) | the points; a named `<Nodes>` block is also a `point` region, as it is a node set in FEBio |
| `<Elements type= name=>` | one cell block, and a `cell` region of that name tagged with its domain's material id (spec 2.5: the `mat` attribute) |
| `<NodeSet>` | a `point` region, in every spec's form: a comma list (4.0 ranges `a:b[:step]` included), `<node id=/>` or `<n id=/>` children, `<node_list>`, or a 2.5 `<NodeSet node_set=/>` include |
| `<ElementSet>` | a `cell` region |
| `<Surface>` | a `side` region when every facet is a face of a solid element, matched by its corner nodes; otherwise its facets become their own cell block and a `cell` region, with a warning |
| `<Edge>` | a `line`/`line3` block and a `cell` region |
| `<DiscreteSet>` | a `line` block (one cell per `<delem>` spring) and a `cell` region |
| `<MeshData>` `<NodeData>` / `<ElementData>` | point / cell data, NaN outside the set they are defined on; the width comes from `data_type` (`scalar`, `vec2`, `vec3`, `mat3s`, `mat3`), or from the values when it is absent (a `type="fiber"` array) |
| `<SurfaceData>`, `<SurfacePair>`, `<PartList>`, materials, loads, boundary conditions, steps | not read |

Element types map to meshio++ cells: `tet4` (and `ut4`), `tet10`, `penta6`, `penta15`, `pyra5`, `pyra13`, `hex8`, `hex20`, `hex27`, `quad4` (and the `q4eas`/`q4ans`/`q4s` shells), `quad8`, `quad9`, `tri3` (and `tri3s`), `tri6`, `tri7`, `tri10`, `line2`/`truss2` and `line3`. Integration-rule aliases such as `TET10G4` or `HEX8G1` name the same nodes. FEBio's `tet5` and `tet15` have no meshio++ cell, so they are a `ReadError`, or with `lenient=True` (`ReadOptions::mLenient`) are read as `tetra`/`tetra10` from their leading nodes, with a warning. `tet20` is always refused.

## Node order

Only `hex27` differs from meshio++'s (VTK's) order. FEBio puts the four mid-height face centres at `s=−1, r=+1, s=+1, r=−1`, while VTK has `r=−1, r=+1, s=−1, s=+1`. The `"febio"` table of the [node-ordering registry](../node_ordering.md) handles it; it comes from the shape functions of `FEHex27` in FEBio's `FECore/FESolidElementShape.cpp`. `tet10`, `penta15`, `pyra13`, `hex20`, `tri6` and `quad8` are in VTK order, which is also what FEBio Studio's VTK exporter assumes.

## Writing

The writer emits spec 4.0:

- `<Module type="solid"/>`, then a `<Material>` section with one **placeholder** material per domain (`isotropic elastic`, `E=1`, `v=0.3`, named after the domain). FEBio refuses a domain whose material is undefined, and meshio++ knows no materials, so replace them. A shell domain gets no `shell_thickness` either; FEBio needs one to run.
- `<Mesh>`: one `<Nodes>` block (ids from 1), then one block per cell block:
  - a 3-D block is an `<Elements>` block with a `<SolidDomain>`;
  - a 2-D block whose cells are all faces of solid cells is a `<Surface>`, otherwise `<Elements>` with a `<ShellDomain>`;
  - a line block in a mesh with 2-D or 3-D cells is an `<Edge>` when its lines lie along cell edges, and a `<DiscreteSet>` (springs) when two-node lines join nodes no cell connects; in a mesh of lines only, it is `<Elements>` with a `<BeamDomain>`.
- A block is named after the `cell` region that covers exactly its cells, else `Part<k>`, `Surface<k>`, `Edge<k>` or `DiscreteSet<k>`.
- Other regions: a `point` region is a `<NodeSet>`, a `cell` region an `<ElementSet>` over the `<Elements>` cells it holds, and a `side` region a `<Surface>` (or an `<Edge>` for the edges of a surface mesh). A region's tag is not kept, since the material ids are renumbered.
- The provenance block is an XML comment inside `<febio_spec>`, so the root tag stays within the bytes content sniffing reads.
- **Dropped, with a warning and a provenance note:** vertex cells and data arrays (MeshData is not written yet; see the [roadmap](../roadmap.md)). Polygons, polyhedra and higher-order Lagrange cells are a `WriteError`.

The written file is meant to be imported into FEBio Studio, or pulled into a model file with `<Mesh from="mesh.feb"/>`, which takes only its `<Mesh>`.

## Validation

These checks were made outside the repository, with FEBio 4.12 built from source:

- The three test fixtures, written by `tools/gen_febio_fixtures.py` in FEBio's numbering, are read by FEBio itself.
- febio-python's spec-3.0/4.0 sample models (MIT) read identically in both engines. Their meshes, written back by this writer and pulled into the original models with `<Mesh from=…>`, give FEBio's original displacements: to 1e-13 for a static surface-load model, 1e-10 for a fibre-reinforced one, and within the solver's tolerance (1e-5 of the peak) for a dynamic model with springs.
