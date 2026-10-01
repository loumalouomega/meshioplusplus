# Interfaces and contact

## Region adjacency

`region_adjacency(mesh, regions=None)` returns a mesh containing the conforming facets shared by different selected Cell regions. Facets are matched by their sorted corner-node ids, so this is exact topology matching rather than geometric proximity; separately numbered coincident surfaces do not match.

When `regions` is omitted or empty, every Cell region is selected in canonical `(kind, name, dim, tag)` order. If fewer than two Cell regions exist and the mesh has at least two cell blocks, each block is treated as a group in block order. Otherwise, pass an ordered sequence of region names or selector dictionaries such as `{"name": "solid", "dim": 3, "tag": 4}`; a selector must resolve to exactly one Cell region.

Each matching facet is emitted once for each selected-region pair. A 3-D source cell contributes faces; a 2-D source cell contributes edges. Facets are identified by corner nodes, including for higher-order cells; polygon and polyhedron faces are supported. For a facet used by multiple cells on either side, the parent cell/facet metadata names the first owner in block-major order and `interface:shared_count` records the number of distinct participating source cells across both sides.

The output keeps the source points, point data, field data and Point regions. Cell data records `interface:region_a`/`interface:region_b` (indices into the selected-group order), `interface:parent_cell_a`/`interface:parent_cell_b`, `interface:parent_facet_a`/`interface:parent_facet_b`, `interface:shared_count` and `interface:measure`. The measure is edge length in 2-D and fan-triangulated area in 3-D. Region Cell groups named `adjacency:<a>:<b>` identify output facets by selected-group index; their entry count and the sum of `interface:measure` give the shared-facet count and area for that pair. Source Cell and Side regions are not copied because their entity indices refer to the input mesh.

```python
interfaces = meshioplusplus.region_adjacency(mesh, ["part_a", "part_b"])
length_or_area = interfaces.cell_data["interface:measure"][0].sum()
```

The Python CLI writes the resulting facet mesh and reports the total facet count and measure:

```sh
meshioplusplus region-adjacency input.vtu interfaces.vtu --regions part_a part_b --json
```

The operation is also available through the C++/C APIs, Fortran, Julia, R, WebAssembly, the native CLI, the pipeline (`RegionAdjacency` with a `Regions` list), and MCP. See [the MCP tool](./mcp.md#mesh-operations) and [the pipeline schema](./pipeline.md).

## Find an interface

`find_interface(mesh, region_a, region_b, *, mesh_b=None, mode="conforming", master="a", gap_tolerance=0.0, angle_tolerance=30.0, overlap_tolerance=0.0, return_report=False)` accepts two Cell-region names, `Region` objects, or selector dictionaries. With `mesh_b=None`, both parts are selected from one mesh; otherwise B is selected from a second mesh. A selector can pin `dim` and `tag` when names are ambiguous, and `block:<index>` selects a cell block.

In `conforming` mode, a boundary facet of A matches a boundary facet of B exactly when their corner-node sets agree; separate meshes must therefore already share node ids. In `proximity` mode, the operation projects each A-facet centroid to B, requires opposing facet normals within `angle_tolerance`, and accepts a facet only when its centroid and every corner are within `gap_tolerance + overlap_tolerance` of B. A zero `gap_tolerance` derives one percent of the mean boundary-edge length. The closest B facet wins, with distance then global cell/facet ids as the deterministic tie-break. This is a proximity query, not a mortar overlap or a boolean intersection.

The returned mesh contains the selected master side's facets and `interface:parent_cell`, `interface:parent_facet`, `interface:partner_cell`, `interface:partner_facet`, `interface:gap`, and `interface:measure` cell data. `interface:gap` is signed along the selected master facet's normal. When `return_report=True`, the function returns `(mesh, report)`; the report includes `num_pairs`, `area` (length for 2-D interfaces), `max_gap`, unmatched boundary-facet counts, and `side_a`/`side_b` Side regions. Each Side region's `(global_cell, local_facet)` entries refer to its own source mesh and can be passed to `split_interface` on that mesh. The facet mesh, report and Python CLI are available through Python, C++, the pipeline, and MCP.

```python
interface, report = meshioplusplus.find_interface(
    mesh, "left_part", "right_part", mode="proximity", return_report=True
)
print(report["num_pairs"], report["area"], report["max_gap"])
```

## Contact pairs

`contact_pairs(slave_mesh, slave_points, master_cells, *, master_mesh=None, tolerance=0.0, require_complete=False)` projects every point in the named Point region to the closest queryable master facet in the selected Cell region. A 3-D volume region contributes its boundary; a 3-D surface-cell region contributes its cells; a 2-D planar mesh contributes boundary edges. A zero tolerance derives one percent of the master's mean facet-edge length. `require_complete=True` raises if any slave point has no facet within tolerance.

The result is a dictionary of arrays: `slave_point`, `master_cell`, `master_facet`, `master_subfacet`, `local_coordinates`, `closest_point`, signed `gap`, unit `normal`, and `unmatched`. Triangle coordinates are three barycentric weights; for a polygon triangulated as a fan, `master_subfacet` is the zero-based fan-triangle ordinal needed to interpret those weights. Edge coordinates use the first two slots as endpoint weights and report subfacet zero. Unmatched rows use `-1` ids and zero-valued geometric fields. Equal-distance facets tie by global cell id then local facet id.

## Split an interface

`split_interface(mesh, side, *, add_cohesive=False, return_report=False)` accepts a Side `Region` or the name/selector of one already on the mesh. It removes selected facets from the node-local cell dual graph, labels the remaining incident-cell fans, and duplicates each point for every additional fan. The component containing the most cells named by the selected Side entries keeps the original point; ties go to the component with the lowest global cell id. Copies are appended in original-point then component order. Point data is copied to each duplicate, Point regions include their duplicate points, and existing cell/field data and regions are retained.

With `add_cohesive=True`, each selected two-owner edge/face adds a zero-thickness `line`, `wedge`, or `hexahedron` from 2-D edges, triangles, or quads. A 2-D cohesive line's connectivity is trace A; its two trace-B point ids are stored in Int64 cell data `cohesive:trace_b` with two components. Other newly added cell-data values are zero initialized. Polyhedron inputs are rejected. When the Side region includes both owners of a facet, only one cohesive cell is inserted.

```python
split_mesh, split_report = meshioplusplus.split_interface(
    mesh, report["side_a"], add_cohesive=True, return_report=True
)
```

![A matched triangular interface and a slave-to-master contact projection](/images/interface_contact.png)

Both CLIs provide `region-adjacency`, `find-interface`, `contact-pairs`, and `split-interface`; `find-interface` can read B from a second file, while `split-interface --side-entries` accepts the JSON `(cell, facet)` rows from its report. The pipeline supports `RegionAdjacency`, `FindInterface`, and `SplitInterface` for one-mesh operations. MCP and the C, Fortran, Julia, R and WebAssembly bindings expose all four operations.
