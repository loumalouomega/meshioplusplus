# Supported Formats

## Format table

Each format name links to a detailed reference page (structure, options, data mapping, and the C++ vs Python behaviour). The **Round trip** column summarises what a write followed by a read keeps (cell types, point/cell/field data, region kinds), as checked by the [format conformance matrix](./conformance.md).

| Format name | Extensions | Read | Write | Extra dependencies | Round trip |
|---|---|---|---|---|---|
| [`abaqus`](./formats/abaqus.md) | `.inp` | ✓ | ✓ | — | [6/8 cells · no data · regions CPS](./conformance.md#abaqus) |
| [`abaqus_fil`](./formats/abaqus_fil.md) | `.fil` (ASCII and binary results) | ✓ | — | — | [read-only](./conformance.md#abaqus-fil) |
| [`ansys`](./formats/ansys.md) | `.msh` | ✓ | ✓ | — | [4/8 cells · no data · regions C](./conformance.md#ansys) |
| [`ansysInp`](./formats/ansysinp.md) | `.cdb`, `.inp` | ✓ | ✓ | — | [8/8 cells · no data · regions CP](./conformance.md#ansysinp) |
| [`ansys_rst`](./formats/ansys_rst.md) | `.rst`, `.rth` | ✓ | — | — | [read-only](./conformance.md#ansys-rst) |
| [`ansys_rst_cyclic`](./formats/ansys_rst.md#cyclic-symmetry) | — (by name: the full rotor of a static cyclic `.rst`) | ✓ | — | — | [read-only](./conformance.md#ansys-rst-cyclic) |
| [`avsucd`](./formats/avsucd.md) | `.avs` | ✓ | ✓ | — | [8/8 cells · data PC](./conformance.md#avsucd) |
| [`cae`](./formats/cae.md) | `.npz` | ✓ | ✓ | — | [0/8 cells · no data](./conformance.md#cae) |
| [`cgns`](./formats/cgns.md) | `.cgns` | ✓ | ✓ | `h5py` | [8/8 cells · data P](./conformance.md#cgns) |
| [`code_aster`](./formats/code_aster.md) | `.mail` | ✓ | ✓ | — | [8/8 cells · no data · regions CP](./conformance.md#code-aster) |
| [`dex`](./formats/dex.md) | `.dex` | ✓ | ✓ | — | [0/8 cells · data P](./conformance.md#dex) |
| [`dolfin-xml`](./formats/dolfin.md) | `.xml` | ✓ | ✓ | — | [1/8 cells · data PC](./conformance.md#dolfin-xml) |
| [`elmer`](./formats/elmer.md) | a directory (`mesh.header`, …) | ✓ | ✓ | — | [7/8 cells · no data · regions C](./conformance.md#elmer) |
| [`ensight`](./formats/ensight.md) | `.case` / `.geo` | ✓ | ✓ | — | [8/8 cells · data PC](./conformance.md#ensight) |
| [`exodus`](./formats/exodus.md) | `.e`, `.exo`, `.ex2` | ✓ | ✓ | `netCDF4` | [8/8 cells · data PC · regions CP](./conformance.md#exodus) |
| [`febio`](./formats/febio.md) | `.feb` | ✓ | ✓ | — | [7/8 cells · data PC · regions CPS](./conformance.md#febio) |
| [`femap`](./formats/femap.md) | `.neu` | ✓ | ✓ | — | [8/8 cells · data PC · regions CP](./conformance.md#femap) |
| [`flac3d`](./formats/flac3d.md) | `.f3grid` | ✓ | ✓ | — | [6/8 cells · no data](./conformance.md#flac3d) |
| [`flux`](./formats/flux.md) | `.pf3` | ✓ | ✓ | — | [8/8 cells · no data](./conformance.md#flux) |
| [`frd`](./formats/frd.md) | `.frd` | ✓ | — | — | [read-only](./conformance.md#frd) |
| [`freefem`](./formats/freefem.md) | `.msh` | ✓ | ✓ | — | [2/8 cells · no data](./conformance.md#freefem) |
| [`gid`](./formats/gid.md) | `.post.msh` / `.post.res`, `.post.bin`, `.post.h5` | ✓ | ✓ | *writing* needs zlib (vendored gidpost); *reading* needs nothing for ascii, zlib for binary, HDF5 for hdf5 | [8/8 cells · data PC](./conformance.md#gid) |
| [`gltf`](./formats/gltf.md) | `.glb`, `.gltf` | — | ✓ | — | [write-only](./conformance.md#gltf) |
| [`gmsh` / `gmsh22`](./formats/gmsh.md) | `.msh` | ✓ | ✓ | — | `gmsh` [fails](./conformance.md#gmsh) / `gmsh22` [8/8 cells · data PCF · regions C](./conformance.md#gmsh22) |
| [`h5m`](./formats/h5m.md) | `.h5m` | ✓ | ✓ | `h5py` | [3/8 cells · data P](./conformance.md#h5m) |
| [`hmf`](./formats/hmf.md) | `.hmf` | ✓ | ✓ | `h5py` | [8/8 cells · data PC](./conformance.md#hmf) |
| [`ip`](./formats/ip.md) | `.ip` | ✓ | ✓ | — | [0/8 cells · data P](./conformance.md#ip) |
| [`libmesh`](./formats/libmesh.md) | `.xda`, `.xdr` (and `.gz`, `.bz2`) | ✓ | ✓ | native gzip reading needs zlib; bzip2 and compressed writing are Python-only | [8/8 cells · no data · regions CPS](./conformance.md#libmesh) |
| [`lsdyna`](./formats/lsdyna.md) | `.k`, `.key`, `.dyn` | ✓ | ✓ | — | [8/8 cells · no data · regions CPS](./conformance.md#lsdyna) |
| [`lsdyna_binout`](./formats/lsdyna_binout.md) | none: the file named `binout` (or `binout0000`...) | ✓ | — | — | [read-only](./conformance.md#lsdyna-binout) |
| [`lsdyna_d3plot`](./formats/lsdyna_d3plot.md) | none: the file named `d3plot` or `d3part` (and its `d3plot01`... family) | ✓ | — | — | [read-only](./conformance.md#lsdyna-d3plot) |
| [`marc`](./formats/marc.md) | `.dat` (a Marc input deck; Tecplot's otherwise) | ✓ | ✓ | — | [7/8 cells · no data · regions CP](./conformance.md#marc) |
| [`marc_t19`](./formats/marc.md#the-post-file) | `.t19` (formatted post file) | ✓ | — | — | [read-only](./conformance.md#marc-t19) |
| [`mdpa`](./formats/mdpa.md) | `.mdpa` | ✓ | ✓ | — | [8/8 cells · data P](./conformance.md#mdpa) |
| [`med`](./formats/med.md) | `.med` | ✓ | ✓ | `h5py` | [8/8 cells · data PC · regions CP](./conformance.md#med) |
| [`medit`](./formats/medit.md) | `.mesh`, `.meshb` | ✓ | ✓ | — | [7/8 cells · no data](./conformance.md#medit) |
| [`mfem`](./formats/mfem.md) | `.mesh` (by content), `.gf` grid functions | ✓ | ✓ | — | [6/8 cells · no data · regions C](./conformance.md#mfem) |
| [`mff`](./formats/mff.md) | `.mff` | ✓ | ✓ | — | [0/8 cells · no data](./conformance.md#mff) |
| [`mfm`](./formats/mfm.md) | `.mfm` | ✓ | ✓ | — | [fails](./conformance.md#mfm) |
| [`mphbin`](./formats/mphbin.md) | `.mphbin` | ✓ | ✓ | — | [8/8 cells · no data · regions C](./conformance.md#mphbin) |
| [`mphtxt`](./formats/mphtxt.md) | `.mphtxt` | ✓ | ✓ | — | [8/8 cells · no data · regions C](./conformance.md#mphtxt) |
| [`nastran`](./formats/nastran.md) | `.bdf`, `.fem`, `.nas` | ✓ | ✓ | — | [8/8 cells · no data · regions C](./conformance.md#nastran) |
| [`nastran_h5`](./formats/nastran_h5.md) | `.h5` | ✓ | — | `h5py` | [read-only](./conformance.md#nastran-h5) |
| [`nastran_op2`](./formats/nastran_op2.md) | `.op2` | ✓ | — | — | [read-only](./conformance.md#nastran-op2) |
| [`netgen`](./formats/netgen.md) | `.vol`, `.vol.gz` | ✓ | ✓ | — | [8/8 cells · no data](./conformance.md#netgen) |
| [`neuroglancer`](./formats/neuroglancer.md) | (no extension) | ✓ | ✓ | — | [1/8 cells · no data](./conformance.md#neuroglancer) |
| [`obj`](./formats/obj.md) | `.obj` | ✓ | ✓ | — | [2/8 cells · no data](./conformance.md#obj) |
| [`off`](./formats/off.md) | `.off` | ✓ | ✓ | — | [2/8 cells · no data](./conformance.md#off) |
| [`openfoam`](./formats/openfoam.md) | `.foam` | ✓ | ✓ | — | [4/8 cells · no data · regions CPS](./conformance.md#openfoam) |
| [`patran`](./formats/patran.md) | `.pat`, `.out` | ✓ | ✓ | — | [7/8 cells · no data · regions CP](./conformance.md#patran) |
| [`pcd`](./formats/pcd.md) | `.pcd` | ✓ | ✓ | — | [1/8 cells · data P](./conformance.md#pcd) |
| [`permas`](./formats/permas.md) | `.post`, `.post.gz`, `.dato`, `.dato.gz` | ✓ | ✓ | — | [8/8 cells · no data](./conformance.md#permas) |
| [`ply`](./formats/ply.md) | `.ply` | ✓ | ✓ | — | [0/8 cells · no data](./conformance.md#ply) |
| [`pmsh`](./formats/pmsh.md) | `.pmsh` | ✓ | ✓ | — | [1/8 cells · data PCF](./conformance.md#pmsh) |
| [`pvd`](./formats/pvd.md) | `.pvd` | ✓ | ✓ | — | [8/8 cells · data PCF](./conformance.md#pvd) |
| [`pvtp`](./formats/pvtp.md) | `.pvtp` | ✓ | ✓ | — | [4/8 cells · data PCF](./conformance.md#pvtp) |
| [`pvtu`](./formats/pvtu.md) | `.pvtu` | ✓ | ✓ | — | [8/8 cells · data PCF](./conformance.md#pvtu) |
| [`radioss`](./formats/radioss.md) | `.rad` (starter deck) | ✓ | ✓ | — | [7/8 cells · no data · regions CPS](./conformance.md#radioss) |
| [`radioss_anim`](./formats/radioss_anim.md) | `<run>A001`… (animation files, by name or content) | ✓ | — | — | [read-only](./conformance.md#radioss-anim) |
| [`radioss_th`](./formats/radioss_th.md) | `<run>T01`… (time-history files, by name or content) | ✓ | — | — | [read-only](./conformance.md#radioss-th) |
| [`stl`](./formats/stl.md) | `.stl` | ✓ | ✓ | — | [0/8 cells · no data](./conformance.md#stl) |
| [`su2`](./formats/su2.md) | `.su2` | ✓ | ✓ | — | [6/8 cells · no data](./conformance.md#su2) |
| [`svg`](./formats/svg.md) | `.svg` | — | ✓ | — | [write-only](./conformance.md#svg) |
| [`szplt`](./formats/szplt.md) | `.szplt` (Tecplot SZL, also by content) | ✓ | — | TecIO (a core built with it, or a shared TecIO) | [read-only](./conformance.md#szplt) |
| [`tecplot`](./formats/tecplot.md) | `.dat` (unless it is a Marc deck), `.tec`, `.plt` (binary, read only) | ✓ | ✓ | — | [5/8 cells · data PC · regions C](./conformance.md#tecplot) |
| [`tetgen`](./formats/tetgen.md) | `.ele` / `.node` | ✓ | ✓ | — | [1/8 cells · no data](./conformance.md#tetgen) |
| [`tikz`](./formats/tikz.md) | `.tikz` | — | ✓ | — | [write-only](./conformance.md#tikz) |
| [`triangle`](./formats/triangle.md) | `.node` / `.ele` / `.poly` | ✓ | ✓ | — | [1/8 cells · no data](./conformance.md#triangle) |
| [`ugrid`](./formats/ugrid.md) | `.ugrid` | ✓ | ✓ | — | [6/8 cells · no data](./conformance.md#ugrid) |
| [`unv`](./formats/unv.md) | `.unv`, `.uff` | ✓ | ✓ | — | [7/8 cells · data PC · regions CP](./conformance.md#unv) |
| [`usd`](./formats/usd.md) | `.usd`, `.usda`, `.usdc` | ✓ | ✓ | `usd-core` | [0/8 cells · no data](./conformance.md#usd) |
| [`vti`](./formats/vti.md) | `.vti` | ✓ | ✓ | — | [fails](./conformance.md#vti) |
| [`vtk` / `vtk42` / `vtk51`](./formats/vtk.md) | `.vtk` | ✓ | ✓ | — | [8/8 cells · data PC](./conformance.md#vtk) |
| [`vtkhdf`](./formats/vtkhdf.md) | `.vtkhdf`, `.hdf` | ✓ | ✓ | `h5py` | [8/8 cells · data PCF](./conformance.md#vtkhdf) |
| [`vts`](./formats/vts.md) | `.vts` | ✓ | ✓ | — | [fails](./conformance.md#vts) |
| [`vtr`](./formats/vtr.md) | `.vtr` | ✓ | ✓ | — | [fails](./conformance.md#vtr) |
| [`vtm`](./formats/vtm.md) | `.vtm` | ✓ | ✓ | — | [8/8 cells · data PC](./conformance.md#vtm) |
| [`vtp`](./formats/vtp.md) | `.vtp` | ✓ | ✓ | — | [4/8 cells · data PCF](./conformance.md#vtp) |
| [`vtu`](./formats/vtu.md) | `.vtu` | ✓ | ✓ | — | [8/8 cells · data PCF](./conformance.md#vtu) |
| [`vtx`](./formats/vtx.md) | `.bp` (DOLFINx VTX: an ADIOS2 directory, also by content) | ✓ | — | ADIOS2 (a core built with it, or `adios2`) | [read-only](./conformance.md#vtx) |
| [`wkt`](./formats/wkt.md) | `.wkt` | ✓ | ✓ | — | [1/8 cells · no data](./conformance.md#wkt) |
| [`xdmf`](./formats/xdmf.md) | `.xdmf`, `.xmf` | ✓ | ✓ | `h5py` (for HDF data) | [8/8 cells · data PC](./conformance.md#xdmf) |
| [`xplt`](./formats/xplt.md) | `.xplt` | ✓ | — | — (zlib for compressed files) | [read-only](./conformance.md#xplt) |
| [`xyz`](./formats/xyz.md) | `.xyz`, `.xyzn`, `.xyzrgb`, `.asc`, `.pts`, `.txt` | ✓ | ✓ | — | [1/8 cells · data P](./conformance.md#xyz) |
| [`z88`](./formats/z88.md) | `z88i1.txt`, `z88structure.txt` (by file name), results `z88o2.txt`/`z88o3.txt` | ✓ | ✓ (structure file) | — | [3/8 cells · no data · regions CP](./conformance.md#z88) |
| [`zarr`](./formats/zarr.md) | `.zarr` | ✓ | ✓ | `zarr` (writing needs 3.x) | [1/8 cells · data PCF](./conformance.md#zarr) |

**Note on directory formats:** [`elmer`](./formats/elmer.md), `openfoam`, [`pmsh`](./formats/pmsh.md) and [`zarr`](./formats/zarr.md) write a *directory* rather than a file. Extension dispatch still works (`case.pmsh` and `case.zarr` carry their suffix on the directory name), but a write target with no extension needs an explicit `file_format=`, and none of the four can be read from or written to a buffer. **Reading sniffs a directory by the files it holds** (v16.2.0): a `mesh.header` (or a `partitioning.N` of `part.n.*` files) makes it `elmer`, and a `constant/polyMesh` or `polyMesh` with `owner` and `faces` (or a decomposed `processor0`, or a multi-region `constant/regionProperties`) makes it `openfoam`, so `read("case")` and `convert case out.vtu` need no format; a directory matching both, or neither, is not guessed. A glob over such a set — `read_sequence("out_*.pmsh")` — matches the suffixed ones, which an ordinary file glob would not; an extension-less Elmer directory has to be listed explicitly.

**Note on the physics-ML formats:** [`pmsh`](./formats/pmsh.md), [`zarr`](./formats/zarr.md), [`cae`](./formats/cae.md) and [`usd`](./formats/usd.md) are **Python-only**. They are not in the shared C++ dispatch registry, so they are absent from the WASM, C, Fortran, Julia, R and native-CLI surfaces; everything else in this table is reachable from all of them. `pmsh`, `zarr` and `cae` are also *lossy by design* — each reduces a mesh to what its consumer's data model holds (one simplex kind, or a triangulated skin plus node fields) — so they are export targets rather than interchange formats.

**Note on LS-DYNA (`lsdyna`)** (v15.2.0): a keyword deck read in full by both engines, following `*INCLUDE` and `*INCLUDE_PATH`. `*PART` becomes a cell region (title as name, `pid` as tag) and `*SET_NODE` / `*SET_SOLID`, `_SHELL`, `_BEAM`, `_PART` / `*SET_SEGMENT` become point, cell and side regions. The standard, `LONG=`, `I10=` and comma-separated card formats are read per card, so they can be mixed in one file, and the tetra, pyramid and wedge that LS-DYNA writes as hexahedra with repeated nodes are collapsed on read and expanded on write. Only geometry is read — materials, sections, contacts and loads are skipped — and the writer puts placeholder section and material ids on every part. `.k`, `.key` and `.dyn` resolve to `lsdyna` (none was claimed before). See [LS-DYNA](./formats/lsdyna.md).

**Note on Elmer meshes (`elmer`)** (v16.2.0): ElmerSolver's native mesh, a *directory* of `mesh.header`, `mesh.nodes`, `mesh.elements`, `mesh.boundary` and `mesh.names`, read and written by both engines. Bulk and boundary elements are separate cell blocks; bodies and boundaries are `cell` regions tagged with their Elmer ids and named from `mesh.names`. Only the `820`/`827` bricks have a node permutation, pinned against ElmerSolver's `elements.def` and its own VTU writer. A partitioned mesh (`partitioning.N/part.n.*`, halo copies included) is merged with each cell's part in `cell_data["partition:part"]`, and `piece=` reads one part. ElmerGrid's binary `mesh.*.bin` files are read too (v16.10.0). The writer writes the highest-dimensional cells as the bulk, every lower-dimensional cell and every `side` region facet as a boundary element with regenerated parents, and, for a mesh with `partition:part`, a `partitioning.N` directory ElmerSolver runs in parallel (v16.10.0), with `ElmerGrid -halo`'s halo elements on request (`halo=True`, v16.11.0). Checked against ElmerGrid and ElmerSolver built from source. See [Elmer](./formats/elmer.md).

**Note on Ansys MAPDL (`ansysInp`, `ansys_rst`)** (v16.3.0): the coded database (`.cdb`) is read and written by both engines with each block cut by its own Fortran format line, so MAPDL, Workbench and HyperMesh decks read alike. Elements become cells by their element routine's category, with degenerate bricks and shells resolved by their repeated nodes and missing midside nodes created at their edge midpoints; the routine, type, material, real constant and section are `ansys:*` cell data, and `CMBLOCK` components are point and cell regions. The writer emits MAPDL's own layouts (wedges, pyramids and tetrahedra as degenerate SOLID185/186 bricks) and keeps a model's element types. MAPDL's binary results (`.rst`, `.rth`) are **read-only**: the same mesh from the geometry records, and one step per result set whose nodal DOF solution (`U`, `ROT`, `TEMP` ...) is point data rotated to the global axes. Since v16.8.0 each set also gives the reaction forces (`RF`, `RMOM`), and the element nodal stresses and strains (`S`, `EPEL` ...) and forces (`ENF`): per element node as cell data, and stresses and strains averaged at the corner nodes as point data. The main file of a distributed solve reads its partial files with it, and `ansys_rst_cyclic` expands a static cyclic-symmetry model to the full rotor. Both are checked against mapdl-archive and pymapdl-reader on files MAPDL wrote. See [`ansysinp.md`](./formats/ansysinp.md) and [`ansys_rst.md`](./formats/ansys_rst.md).

**Note on FEBio (`febio`, `xplt`)** (v16.2.0): the mesh of an FEBio input file (`.feb`, spec 2.5 in `<Geometry>`, 3.0 and 4.0 in `<Mesh>`) is read by both engines under the rules of FEBio's own parsers: each `<Elements>` block a cell block and a cell region, `<NodeSet>`/`<ElementSet>` point and cell regions, a `<Surface>` on solid faces a side region, `<Edge>`/`<DiscreteSet>` line blocks, `<MeshData>` point and cell data. The writer emits spec 4.0 with a placeholder material per domain; faces of solids become surfaces, lone two-node lines discrete springs, point and cell data `<MeshData>` on `meshdata:<name>` sets (v16.10.0). FEBio's plot file (`.xplt`, read-only) is read one state at a time (`time_step`), compressed or not, in either byte order: nodal variables as point data, per-element ones as cell data, per-element-node ones averaged to the points; since v16.10.0 each surface is a facet block carrying its surface variables, and a remeshed run reads each state on its own mesh; since v16.11.0 each edge is a line block carrying its edge variables. Only `hex27` has a node permutation (`FEHex27`'s face centres). Both were checked against FEBio 4.12 built from source. See [FEBio input](./formats/febio.md) and [FEBio plot files](./formats/xplt.md).

**Note on Code_Aster meshes (`code_aster`)** (v16.0.0): Code_Aster's own ASCII mesh, read and written by both engines under the rules of Code_Aster's reader: only the first 80 columns of a line are read, records are token streams that may wrap, `%` starts a comment. `COOR_nD` gives the points, `POI1` … `HEXA27` the cell blocks (`TRIA7` as the new `triangle7`), and `GROUP_MA`/`GROUP_NO` cell and point regions without a tag. Its node order is **not** MED's: the quadratic hexahedra and wedges list the vertical mid-edges before the top ring, and the tables are pinned against Code_Aster's own gmsh and MED readers in the [node-ordering registry](./node_ordering.md). The writer keeps every line within 80 columns, names nodes and elements `N…`/`M…`, sanitises group names to 24 characters, and drops side regions and data arrays with a warning. See [Code_Aster](./formats/code_aster.md).

**Note on the Patran and Femap neutral files (`patran`, `femap`)** (v16.5.0): two FEM interchange files, read and written by both engines. Patran's `.pat`/`.out` is a sequence of fixed-width `(I2,8I8)` packets: nodes, elements (shape in the header, linear or quadratic by node count, `hexahedron20`/`wedge15` with the vertical mid-edges before the top ring), and named components as point and cell regions; elements no component names are grouped by property. Femap's `.neu` is a sequence of `-1`-delimited, comma-separated blocks whose layouts change with the version, so every record is read by position and length (4.41 to 2020.1): nodes, elements in Femap's 20-slot degenerate-brick node layout, property titles, groups as regions, and output sets as steps (`time_step`) whose `451`/`1051` vectors become point or cell data. The Femap writer emits the 8.2 layout: the mesh, groups and, since v16.11.0, results as an output set. Since v16.12.0 Patran's loads and boundary conditions (packets 06, 07, 08, 10 and 11) are read and written as `patran:` data, and Patran 2.5 result files (`.nod`/`.dis`/`.els`, text or binary) are read onto the mesh with `patran.read(path, results=...)`. Neither Patran nor Femap was available: the Patran fixtures are written from its documentation, the Femap ones are real Femap, EMSolution and MYSTRAN files, checked against FrontISTR's `neu2fstr` and femap_neutral_parser. See [Patran](./formats/patran.md) and [Femap](./formats/femap.md).

**Note on MFEM meshes (`mfem`)** (v16.5.0): MFEM's own `.mesh` (v1.0–v1.3, conforming) and its `.gf` grid functions, read and written by both engines. `.mesh` stays Medit's extension; a file whose first line names an MFEM mesh is read as one, by the Python reader loop and by the native resolver alike. Attributes become `mfem:attribute` and `attribute_<n>`/`boundary_<n>` regions, v1.3 attribute sets named regions. Order-2 `H1` nodes give quadratic cells whose points are MFEM's degrees of freedom in MFEM's own numbering; higher orders keep the vertices, with a warning. Grid functions (`mfem.read(path, {name: gf})`, `mfem.write(..., grid_functions=True)`) become point or cell data. Checked against MFEM 4.10 (PyMFEM) both ways. See [MFEM](./formats/mfem.md).

**Note on libMesh, Z88, Abaqus results and Radioss (`libmesh`, `z88`, `abaqus_fil`, `radioss`)** (v16.7.0): four FEM formats read by both engines. libMesh's `.xda`/`.xdr` (every version from 0.7.0 to 1.8.0, ASCII and XDR) gives the active cells of a refined mesh, subdomains, side sets (carried to refined children) and node sets. Z88's structure file is recognised by its fixed name (`z88i1.txt` before `.txt`'s xyz), read with its `z88o2.txt` displacements and `z88o3.txt` stresses, and written back in the Z88OS v15 layout; Z88R solved the written decks. The Abaqus `.fil` results file (binary in either byte order, or ASCII) gives the model, node and element sets, and one step per increment with nodal and element results; checked against pybaqus on real Abaqus 2023 output. The OpenRadioss starter deck (`*_0000.rad`) gives parts, subsets, groups and surfaces as regions, degenerate bricks as the solids they stand for, and follows `#include`; checked on OpenRadioss's 81 QA decks. Since v16.17.0 it is written too, by both engines alike, and OpenRadioss's own starter accepts the decks written with placeholder materials and properties (`stubs=True`). See [libMesh](./formats/libmesh.md), [Z88](./formats/z88.md), [Abaqus `.fil`](./formats/abaqus_fil.md) and [Radioss](./formats/radioss.md).

**Note on the v16.11.0 completions (`tecplot`, `elmer`, `xplt`, `mfem`, `femap`, `patran`, `libmesh`, `z88`, `radioss`, `radioss_anim`)**: the halves those formats still lacked. Tecplot reads and writes face-based `FEPOLYGON`/`FEPOLYHEDRON` zones (polygon and polyhedron cells). The Elmer writer adds halo elements (`halo=True`, as `ElmerGrid -halo`). FEBio `.xplt` reads edge sections and their variables. MFEM reads and writes meshes and fields of any order as VTK Lagrange cells, reads non-conforming meshes as their leaves and parallel runs as one mesh. The Femap writer writes output sets and the reader takes Femap 2401's property records; Patran reads and writes 9-node quadrilaterals and 7-node triangles, checked on real P3/PATRAN and CUBIT exports. libMesh is written (1.8.0 layout) and read gzip- or bzip2-compressed, with edge and shell-face sets. Z88 reads and writes the whole deck (constraints, materials, element parameters, integration orders) and reads every element family's stresses and the nodal forces, plus Z88Aurora's sets. Radioss applies units and resolves boxes, generators and every surface form, and `radioss_anim` reads the animation files as a transient sequence. See each format's page.

**Note on the v16.12.0 completions (`tecplot`, `mfem`, `patran`, `abaqus_fil`, `radioss`, `radioss_th`, `z88`, `libmesh`, `marc`, `marc_t19`, `ansys_rst`, `lsdyna_d3plot`, `lsdyna_binout`, `nastran_op2`)**: the last halves of the FEM interchange formats, each checked against the files its tool writes where one could be run or found. Tecplot reads binary `.plt` back to version 7 (`#!TDV71`); MFEM reads NURBS meshes, Bernstein and serendipity spaces, grid functions on non-conforming meshes and parallel non-conforming runs; Patran reads and writes loads and boundary conditions (packets 06–08, 10, 11) and reads its `.nod`/`.els`/`.dis` result files; the Abaqus `.fil` reader takes modes, energies, contact output, element matrices and rebar; Radioss applies per-card units, skewed boxes and analytical surfaces, reads 4.x decks and engine decks, and a new `radioss_th` reads its time-history files; Z88 reads the 16-node plate whole and its surface loads, and writes Z88Aurora's sets; libMesh writes refinement trees back and compresses natively (bzip2 with `MESHIOPLUSPLUS_WITH_BZIP2`); Marc follows `INCLUDE`, reads Herrmann, rebar and composite elements, edge and face sets and remeshing increments; the `.rst` reader adds element energies, fluxes, gradients, nonlinear and contact data and the modal cyclic expansion; d3plot reads SPH, airbag particles, rigid bodies and roads, 20- and 27-node solids and `d3part`, and a new `lsdyna_binout` reads the LSDA binout; OP2 reads complex, random and SORT2 tables, element forces, energies, springs and dampers. What still needs a vendor licence to check is listed in the [roadmap](./roadmap.md#awaiting-a-licensed-run).

**Note on DOLFINx VTX and Tecplot SZL (`vtx`, `szplt`)** (v16.13.0): two read-only formats whose only reader is a heavy or vendor library, so neither is ever a default. `vtx` reads DOLFINx's `VTXWriter` output (an ADIOS2 `.bp` directory) as a time sequence, through a core built with `MESHIOPLUSPLUS_WITH_ADIOS2` or the `adios2` package, and matches ParaView's VTX reader on real DOLFINx runs; `szplt` reads Tecplot's subzone-loadable files through TecIO (a core built with `MESHIOPLUSPLUS_WITH_TECIO`, or a shared TecIO through ctypes) exactly as the same data in `.plt`. The formats that only their vendor's software reads — Abaqus `.odb`, Marc `.t16`, Femap `.modfem`, results beyond the `.rst` reader — have [vendor-runtime routes](./vendor_routes.md) instead: scripts that run in the vendor's Python and write files meshio++ reads.

**Note on CalculiX results (`frd`)** (v15.3.0, binary layout and `.dat` print v15.5.0): the result file of `ccx`, **read-only** — it is what a solver writes, so there is nothing to write back. Every `100C` increment is a step for the [sequence engine](./sequences.md) (`time_step`, `read_metadata(...).time_values`, `read_sequence`), each result block is a point data array named as the file names it, and the six components of a symmetric tensor keep the file's order `xx yy zz xy yz zx`. `ccx` expands shells and beams into solids, so the mesh read is not the `.inp` mesh; the he20, pe15 and be3 node orders are permuted to the Abaqus order. The short (`I5`) and long (`I10`) ASCII layouts and the binary layout (`*NODE OUTPUT`/`*ELEMENT OUTPUT`) are all read. `derived=True` on `meshioplusplus.frd.read` adds von Mises and principal values beside each stress and strain tensor, now backed by the shared [`tensor_invariants`](./tensor_invariants.md) operation. `meshioplusplus.frd.read_dat` separately parses the companion `.dat` tabular print (`*NODE PRINT`/`*EL PRINT`) into plain tables, Python-only. See [CalculiX results](./formats/frd.md).

**Note on MSC Nastran HDF5 results (`nastran_h5`)** (v15.7.0): the result database MSC Nastran writes with `MDLPRM,HDF5`, **read-only**. The model comes from `/NASTRAN/INPUT` (GRID points moved to the basic system since v16.10.0, one cell block per element card and order, `nastran:eid`/`nastran:pid` cell data, one region per property id), and every result domain (subcase, mode, time or frequency) that an INDEX table references is a step for the [sequence engine](./sequences.md). Nodal tables become point data (`DISPLACEMENT`, `DISPLACEMENT_ROT`, `EIGENVECTOR_real`...), element tables cell data named `<group>:<member>` (`STRESS:X`), with the corner, ply and station values and the grid point forces as `(cells, columns)` arrays beside them (v16.10.0). Other vendors' HDF5 schemas are refused rather than misread. Checked bit for bit against pyNastran on real MSC Nastran 2020 output. See [MSC Nastran HDF5](./formats/nastran_h5.md).

**Note on LS-DYNA and Nastran binary results (`lsdyna_d3plot`, `lsdyna_binout`, `nastran_op2`)** (v16.9.0): two solver result files, **read-only**, by both engines, every state or (subcase, mode, time) a step for the [sequence engine](./sequences.md). The d3plot family (`d3plot`, `d3plot01`...) is found by its file name or its control block, in single or double precision and either byte order; parts become regions, the degenerate solids collapse as in the keyword reader, per-integration-point and per-layer variables use the `(cells, points * components)` layout, and the deletion flags are the int8 mask `lsdyna:alive`; checked against lasso-python on real and lasso-written families; since v16.12.0 also `d3part`, 20/27-node hexahedra, SPH and airbag particles, rigid roads and rigid-body motion, and the new `lsdyna_binout` reads the LSDA binary output (`nodout` outputs as steps, the other databases as field data). OP2 (MSC and NX, 32- and 64-bit) builds its mesh as `nastran_h5` does, or from the input deck beside a file without geometry, and reads the real SORT1 nodal tables, the stress and strain of rods, bars, beams, shells, composites and solids and the grid point forces under `nastran_h5`'s names and array shapes (coordinate systems, corner, ply and station values since v16.10.0); since v16.12.0 also complex, random and SORT2 tables, element forces, heat fluxes and energies, springs and dampers as cells, `PARAM,POST,-2` files and points from the basic grid point table; checked against pyNastran. See [LS-DYNA d3plot](./formats/lsdyna_d3plot.md) and [Nastran OP2](./formats/nastran_op2.md).

**Note on MSC Marc (`marc`, `marc_t19`)** (v16.8.0): Marc's input deck and formatted post file, read by both engines; the deck is written too since v16.17.0, by name (`file_format="marc"`, since a `.dat` alone writes Tecplot), and the post file stays [read-only by design](./conformance.md#marc-t19). A deck's `COORDINATES` and `CONNECTIVITY` (fixed, `EXTENDED` or free format; Mentat's letterless exponents and touching fields) give the mesh, every Marc element type read being in meshio++'s node order already, and `DEFINE ELEMENT/NODE SET` give cell and point regions (ranges, set names, `AND`/`EXCEPT`/`INTERSECT`). `.dat` stays Tecplot's unless the file opens as a Marc deck, a content check both engines apply. A `.t19` gives the same mesh and sets, and one step per increment: nodal vectors as point data, element post codes per integration point as cell data, a tensor's six codes as one `xx yy zz xy yz zx` array. Checked on real Marc Mentat decks and a real `.t19`; no Marc licence was available. See [Marc](./formats/marc.md).

**Note on the point-cloud formats (`pcd`, `xyz`)** (v15.1.0): both map to the points plus one `vertex` block, exactly what [`subsample_points`](./point_budgets.md) emits, so `.pcd` → `subsample` → `proximity-graph` → `.vtu` runs from the CLI with no intermediate format. `pcd` is PCL's v0.7 layout in all three `DATA` modes (`ascii`, `binary`, `binary_compressed`), with `rgb`/`rgba` unpacked by bit-cast and `normal_x/y/z` gathered into `normals`; `xyz` is a headerless convention, so its columns are resolved from a `columns=` list, a header comment, the column count and extension, or the value ranges, and an ambiguous file is an error rather than a guess. Chemistry XYZ (atom count, comment, `element x y z`) shares the extension and is refused by name. `.txt` and `.asc` now resolve to `xyz`. LAS/LAZ and E57 are out of scope.

**Note on `.msh`:** `ansys`, `freefem`, and `gmsh` all use `.msh`. When writing without an explicit `file_format`, meshio++ picks `gmsh` if the mesh carries gmsh-native tags (`gmsh:physical`/`gmsh:geometrical`/`gmsh:dim_tags`) or MED-derived tags (`cell_tags`/`point_tags`/`med:*`), else falls back to the first registered candidate (`ansys`). When reading, meshio++ tries the registered formats in order and uses the first that parses the file. Specify `file_format` explicitly (e.g. `file_format="freefem"`) to avoid ambiguity either way.

**Note on `.post.*`:** extension resolution tries the *longest* matching suffix first, so a compound extension always wins over a shorter one that also matches — `.post.msh` resolves to `gid`, never falling through to `.msh`'s own `ansys`/`gmsh`/`freefem` candidates (even for a mesh carrying gmsh-native tags, which would otherwise steer a plain `.msh` write to `gmsh`); `.post` alone (no `.msh` suffix) still resolves to `permas`, unaffected.

**Note on `.inp`:** `abaqus` and `ansysInp` both use `.inp`. `abaqus` is registered first, so plain extension-based dispatch resolves to Abaqus by default; pass `file_format="ansysInp"` (or call `meshioplusplus.ansysInp.read`/`write` directly) to select the Ansys/APDL reader for a `.inp` file.

**Note on `tetgen`:** The format spans two files (`.node` + `.ele`). It cannot be read from or written to a buffer.

**Note on `triangle` vs `tetgen` (`.node`/`.ele`):** both formats use `.node`/`.ele`; `tetgen` is registered first, so plain extension dispatch resolves to it. When *reading*, a 2D pair makes tetgen raise and the dispatcher falls through to `triangle` automatically; when *writing*, pass `file_format="triangle"` (or use a `.poly` path, which defaults to `triangle`). Like tetgen, the format spans multiple files and cannot use buffers.

**Note on `ensight`:** EnSight Gold (`.case` + `.geo` sibling pair, ASCII and C-binary with byte-order auto-detection). Multi-part files concatenate into one point array with the owning part recorded as `cell_data["ensight:part"]`. Since v11.3.0 a `.case` file's `TIME`/`VARIABLE` sections are read (`point_data`/`cell_data`/tensors, one step selected by `time_step`); since v15.5.0 they are written too (scalar/vector/tensor symm/tensor asym plus `field_data` as `constant per case`), single-part and non-transient only. Cannot use buffers. The `.geo` extension is also used by Gmsh *script* files, which meshio++ never claimed.

**Note on `vti`:** VTK XML ImageData is a **regular lattice**: its geometry is the `Origin`/`Spacing`/`WholeExtent` attributes rather than a point array. Reading expands the extent into explicit `hexahedron` cells; writing therefore *requires* a lattice, and a mesh that is not one — including a **partial** grid (`voxelize`'s `surface`/`inside` fills, or `compute_sdf`'s octree, whose holes ImageData cannot express) — raises `WriteError` by name. It is the only format that round-trips a generated grid's geometry, which is why [`compute_sdf`](./sdf.md) points at it.

**Note on `vtkhdf`** (v14.0.0): Kitware's HDF5-based VTK format — the one ParaView-native format that keeps time, partitions and fields in a single file. It reads and writes `UnstructuredGrid`, `PolyData`, `PartitionedDataSetCollection` and `MultiBlockDataSet`, including polyhedra, partitioned files (merged into one mesh with one cell region per piece, or `piece=k` for one) and transient files (`time_step=k`; the [time-series writers](./vtkhdf_time_series.md) and the [sequence engine](./sequences.md) write them). The version written is the oldest that covers the content. `.hdf` is registered but is a generic HDF5 extension, so a file without a `/VTKHDF` group is refused with a message asking for an explicit format. See [its own page](./formats/vtkhdf.md) for the layout facts, the composite constraints and the transient-composite limitation.

**Note on `vts`** (v11.6.0): VTK XML StructuredGrid is the same `nx*ny*nz` hexahedron topology as `vti`, except with an explicit `<Points>` array instead of `Origin`/`Spacing` — so, unlike `vti`, **reading never requires a lattice** (a genuinely curved structured mesh reads correctly); writing still does, for the same reason `vti`'s writer does. Degenerate (2-D/1-D) extents are not expanded to quad/line/vertex cells, a documented remainder.

**Note on `vtr`** (v11.6.0): VTK XML RectilinearGrid states the same topology via three independent per-axis coordinate arrays (a tensor product) instead of explicit points. **Reading is fully general** — a genuinely graded (non-uniformly spaced) grid reads correctly, with no uniformity check at all. **Writing still requires a *uniform* lattice**, since recovering three arbitrary per-axis arrays needs a detector that does not exist yet; a graded mesh raises `WriteError` — a documented follow-up.

**Note on `vtm`** (v11.6.0): VTK XML MultiBlock is an index file plus one `.vtu` piece per cell block, so unlike `vti`/`vts`/`vtr` it has **no lattice restriction at all** — any mesh with one or more cell blocks round-trips. Writing carves the mesh into pieces (points pruned per piece via `clean`); reading combines the pieces with [`merge()`](./merge.md) (no welding), attaching one named "cell" region per piece. See [its own page](./formats/vtm.md) for the directory layout and the region-naming asymmetry between read and write.

**Note on `pvd`, `pvtu` and `pvtp`** (v15.0.0): the ParaView index formats, the on-disk face of [`partition`](./partition.md) and of the [sequence engine](./sequences.md). A `.pvtu` (or `.pvtp`, over `.vtp` pieces) declares the arrays every piece holds and names one `.vtu` piece per *part*, where `.vtm` carves by cell *block*: writing carves a mesh by the integer `partition:part` cell array (what `partition_labels` produces), or takes already-carved pieces (`pvtu.write_pieces`, e.g. `partition` output with its `ghost_layers` halos, which become `vtkGhostType` cells and points and the index's `GhostLevel`); every piece must declare identical arrays, checked before anything is written. Reading merges the pieces without welding, one `piece_<i>` cell region each, or keeps one with `piece=`; ghost cells are kept unless `ghosts="drop"` (`--drop-ghosts`, `ReadOptions::mGhosts`, and the same switch on every binding). A `.pvd` is a time-indexed collection of serial or parallel XML files: `timestep=` selects the step (`time_step=`) and `part=` the piece within it (`piece=`), it is a series reader and writer for `read_sequence`/`write_sequence`, and a step that is one `.pvtu` nests without special cases. `.vtu` and `.vtp` carry `field_data` in a `<FieldData>` element (v15.0.0), which is what lets a piece carry a `TimeValue`; a `.pvd` entry with no `timestep=` takes its time from it. ParaView opens all three formats (verified against 6.1.1).

**Note on `vtp`:** VTK XML PolyData holds surface cells only (`vertex`/`line`/`triangle`/`quad`/`polygon`); volume or quadratic cells raise `WriteError`. PolyData has no cell-type array, so 3-/4-noded `polygon` cells read back as `triangle`/`quad`. Triangle strips are not supported.

**Note on `stl` / `ply` and volume meshes:** both writers extract and write the boundary **skin** of a 3D volume mesh by default (`skin=True` — see [Skin extraction](./extract_skin.md); STL additionally triangulates quads, PLY compacts the vertex table). Pass `skin=False` for the legacy behavior (volume cells dropped with a warning).

**Note on `gltf`** (v15.4.0): Write-only. The **surface** of the mesh as glTF 2.0 for the web, three.js and Blender — the skin of volume cells, 2-D cells, `line` cells as `LINES` and `vertex` cells or a cell-less point cloud as `POINTS`, one named node per cell region. `.glb` writes the binary container, `.gltf` JSON with a `.bin` beside it. Normals are per vertex, so points are duplicated at creases (`split_angle`, sharing the smooth-fan kernel of [`compute_normals`](./normals.md)); every one-to-four component `point_data` array is exported raw as `_NAME`; `color_by` bakes a field through a colormap into linear `COLOR_0` with an unlit material. Y-up, `float32`; the axis change, unit scale and recentring offset live on the root node, not in the coordinates. The C++ core and the Python reference write identical bytes and the output passes the Khronos glTF-Validator with zero errors and warnings; see [glTF](./formats/gltf.md).

**Note on `svg`:** Write-only. Flat 2D meshes draw directly; genuinely 3D meshes render their boundary skin (see [Skin extraction](./extract_skin.md)) through an orthographic camera (`azimuth`/`elevation`/`roll`, default the CAD isometric view) with painter's-algorithm depth ordering. C++ core with a Python fallback.

**Note on `tikz`:** Write-only; emits a standalone (directly `pdflatex`-compilable) LaTeX/TikZ document by default (`standalone=False` for a bare `tikzpicture` snippet). Flat 2D meshes draw directly; genuinely 3D meshes render their boundary skin like the SVG writer (same camera parameters). C++ core (byte-identical to the Python reference, including the 3D path) with a Python fallback.

**Note on `openfoam`:** A directory-based format (`points`/`faces`/`owner`/`neighbour`/`boundary` under `constant/polyMesh`), not a single file — so it is the only writer that *creates a directory*. Writing takes a `.foam` marker file, a `polyMesh` directory, or a case root; a case *directory* has no extension, so that form needs an explicit `file_format="openfoam"`. ASCII by default; `binary=True` (v15.5.0) writes little-endian binary at a chosen `label_bits`/`scalar_bits` width, Python/C++ only. Polyhedral cells are native here.

**Note on `mfm`:** Single element type per file (non-hybrid), linear elements only.

**Note on `exodus`:** One-node `SPHERE`/particle meshes (peridynamics solvers such as [PeriLab](https://github.com/PeriHub/PeriLab.jl) write these) read as `vertex` cells, and per-element **attributes** — `attrib{k}`, where a sphere's radius or a shell's thickness lives — round-trip as `cell_data` under the `exodus:attr:` prefix, always float64 and NaN for a block the file gives no such attribute. Since v9.9.0 ordinary (non-attribute) `cell_data` round-trips too, as element variables (`name_elem_var`/`elem_var_tab`/`vals_elem_var{j}eb{k}`); element-block names round-trip through `Cell` regions (`eb_names`), and `field_data["exodus:time"]` labels the written step. See [`exodus.md`](./formats/exodus.md#element-attributes).

**Note on FEconv-derived formats (`unv`, `mfm`, `freefem`, `flux`, `mff`, `dex`, `ip`):** These readers/writers were implemented against the [FEconv](https://github.com/victorsndvg/FEconv) format documentation and public format specs (FEconv is GPL; no FEconv code or data is used — MIT-clean, with fixtures generated by round-trip). `unv` is described on [its own page](./formats/unv.md): since v15.6.0 its node orders and group layouts are checked against gmsh, Salome and pyuff, its groups are regions and its results (2414, 55, 56 and the 58/58b functions of test-lab `.uff` files) are the steps of a sequence; `flux` round-trips per-element region references as `cell_data` (`pf3:ref`), and since v16.6.0 applies FLUX's mirrored solid node order through the [node-ordering registry](./node_ordering.md). `mphtxt` started here too; since v16.1.0 it follows COMSOL's own documentation and files instead, with its binary twin `mphbin` (see [its page](./formats/mphtxt.md)). Node orderings for higher-order elements round-trip losslessly but may differ from the originating tool's internal ordering for some element types. FEconv's own `examples/` (many small meshes, several the same mesh written by different tools) are a useful cross-check but GPL, so they are not committed: `tests/python/test_feconv_examples.py` reads a local checkout when `MESHIOPLUSPLUS_FECONV_DIR` points at it, checking engine parity, orientation, mid-edge placement and the twins (FLUX, VTU, Medit, Gmsh, Fluent, COMSOL and UNV files). The quirks they exposed in v16.6.0 are reproduced by the generated fixtures of `tools/gen_feconv_quirk_fixtures.py`.

**Note on the field-only formats (`mff`, `dex`, `ip`):** These carry result data, not geometry. They read into a geometry-less `Mesh` (no cells) with the field(s) in `point_data`; `dex`/`ip` also populate `points` from the coordinates in the file, while `mff` carries no coordinates (its `points` has zero columns and only the field values round-trip). To attach a field to a mesh, read the field file and the mesh file separately and copy the field `Mesh`'s `point_data` onto the geometry `Mesh` — there is no fixed naming convention pairing a field file with its mesh (unlike TetGen's `.node`/`.ele`).

---

## Provenance

Since v10.15.0, every writer that carries a free-text header slot emits one canonical line, `Written by meshio++ v<release>` (`meshioplusplus._provenance.TAG` on the Python side, `detail::kProvenanceTag` on the C++ side — a single owner on each side, so the two engines emit character-identical bytes). Since v10.16.0 a caller can additionally opt into a richer block — source, target, the operation chain, conversion assumptions, a timestamp — rendered where a format's header slot can hold it; see [`doc/provenance.md`](./provenance.md) for the design and [`doc/roadmap.md`](./roadmap.md) section 1 for what remains open (reading the block back on read).

Before this, roughly two dozen writers hand-wrote their own version of the line and had drifted three ways: a stale `meshio` (not `meshio++`) name in a handful of Python writers, a `(C++ core)`-vs-`v{version}` split that made the C++/Python fallback boundary visible in the output bytes, and three writers (`obj`, `ply`, `exodus`) embedding a wall-clock timestamp that made writing the same mesh twice produce different bytes.

The table below is the audit this fixed: every format's comment syntax (if any), where a header may legally sit, and whether meshio++'s writer uses that slot today. A "—" in the last column means the format admits no free-text slot at all, or the slot it does admit (an ID/label field, a fixed keyword line) is a structural record rather than a place for a human-readable credit, and meshio++ does not write one there. The same classification, at the finer `SlotTier` granularity the opt-in block renders against (`None`/`Bounded`/`SingleLine`/`Block`), is [`doc/provenance.md`](./provenance.md#slot-tiers-and-the-moderequired-interaction)'s table -- this one stays the format-level audit.

| Format | Comment syntax admitted | Where a header may sit | Emits the tag? |
|---|---|---|---|
| `abaqus` | `**`-prefixed lines; free text inside `*HEADING`'s data lines | Anywhere; writer uses `*HEADING` | Yes |
| `ansys` | Section `(0 "text")` is the format's own comment section | Anywhere before `(0 ...)` sections | Yes (writer uses the `(1 "...")` program/version record, not a `(0 ...)` comment) |
| `ansysInp` | `!` or `/` prefix | Anywhere | — |
| `avsucd` | `#` prefix | Top of file | Yes |
| `cae` | A bytes array written as the archive's **first** zip member; NUL-padded rows | The file's first bytes, so the ordinary head scanner finds it | Yes |
| `cgns` | None (HDF5/CGNS binary container) — an HDF5 root attribute is the nearest equivalent | n/a | — |
| `code_aster` | `%` prefix | Anywhere; the writer puts it first | Yes |
| `dex` | None — a fixed two-line header, the second ending in `#` | n/a (structural, not free text) | — |
| `dolfin-xml` | XML `<!-- -->` | Anywhere in the document | — |
| `elmer` | `!` lines in `mesh.names` (any line without both `$` and `=` is ignored by ElmerGrid and ElmerSolver) | `mesh.names`, which the writer always writes; the writer puts it first | Yes |
| `ensight` | None named, but the `.geo` header's description line 2 is free text | Fixed line 2 of the `.geo` header | Yes |
| `exodus` | netCDF root `title` attribute | n/a (one attribute, not a line) | Yes |
| `femap` | None (a block's title line is the nearest free-text slot) | Block 100's title line | Yes (the title) |
| `febio` | XML `<!-- -->` | Anywhere in the document; the writer puts it just inside `<febio_spec>`, so the root tag stays within the bytes sniffing reads | Yes |
| `flac3d` | `*` prefix | Top of file | Yes |
| `flux` | None named — an unlabeled free-text line the reader skips (keyed lookup, not positional) | Top of file | Yes |
| `freefem` | None (fixed positional numeric header) | n/a | — |
| `gid` | `# Name: value` "user attribute" lines (ascii/binary); an HDF5 group attribute (`hdf5` flavour) | One per mesh/result "block", via `GiD_fWriteMeshUserAttribute`/`GiD_fWriteResultUserAttribute` | Yes (ascii confirmed; the `hdf5` flavour's result-side attribute placement is not independently verified — see `gid.md`) |
| `gltf` | `asset.generator` (the tag) and, with a scope open, `asset.extras["meshioplusplus:provenance"]` (the full block) | The `asset` object | — |
| `gmsh` / `gmsh22` | `$Comments`/`$EndComments` section (spec-legal; our reader skips it, neither writer emits one) | Anywhere between sections | — |
| `h5m` | None (HDF5 container) — an HDF5 root attribute is the nearest equivalent | n/a | — |
| `hmf` | None (HDF5 container) — an HDF5 root attribute is the nearest equivalent | n/a | — |
| `ip` | None (fixed positional numeric layout) | n/a | — |
| `lsdyna` | `$` prefix | Anywhere; the writer puts it after `*KEYWORD` | Yes |
| `marc` | `$` prefix | Anywhere; the writer puts it after `TITLE` (v16.17.0) | Yes |
| `mdpa` | `//` prefix (Kratos/C++ style) | Anywhere | — |
| `med` | HDF5 `DES` mesh-description field (user-overridable via `MedInfo.description`, not a provenance slot) | n/a | — (see the `med.md` note below) |
| `medit` | `#` prefix | Anywhere | — |
| `mfem` | `#` lines (between sections) | Anywhere between sections; the writer puts it after the header line | Yes |
| `mff` | None (fixed positional numeric layout, no geometry) | n/a | — |
| `mfm` | None (fixed positional numeric layout) | n/a | — |
| `mphbin` | None (binary serialisation, no comment slot) | n/a | — |
| `mphtxt` | `#` prefix | Top of file | Yes |
| `nastran` | `$` prefix | Top of file, before `BEGIN BULK` | Yes (see the sentinel note below) |
| `netgen` | `#` prefix | Top of file | Yes |
| `neuroglancer` | None (binary chunked octree format) | n/a | — |
| `obj` | `#` prefix | Top of file | Yes |
| `off` | `#` prefix | After the `OFF` magic line | Yes |
| `openfoam` | C-style `/* ... */`; the `FoamFile` banner's fixed-width credit cell is this writer's own convention | Top of file, inside the banner box | Yes (C++ writer only — no Python twin) |
| `patran` | None (the `25` title packet's 80-column card is the nearest free-text slot) | The title card | Yes (the title) |
| `pcd` | `#` prefix (the header's own comment syntax) | After the `# .PCD v0.7` first line, before `VERSION` | Yes |
| `permas` | `!` prefix | Top of file | Yes |
| `ply` | `comment ` prefix (the format's own keyword) | Anywhere in the header, before `end_header` | Yes |
| `pmsh` | None — a memmap directory has no free-text field, and a sidecar would add a file upstream's own loader never wrote | n/a | — |
| `stl` | Binary: an 80-byte free header slot. ASCII: none (the `solid` line names the object, not a comment) | Binary: the first 80 bytes. ASCII: none | Yes (binary only) |
| `su2` | `%` prefix | Anywhere | — |
| `usd` | The root layer's `documentation` metadata field | Stage metadata, at the top of a `.usda` | Yes (a binary crate layer needs `usd-core` to read it back) |
| `svg` | XML `<!-- -->` | Anywhere in the document | — |
| `tecplot` | None named — the `TITLE = "..."` record is the nearest free-text slot | Top of file | Yes |
| `tetgen` | `#` prefix | Top of each file (`.node`/`.ele`/`.poly`) | Yes |
| `tikz` | `%` prefix (LaTeX/TikZ convention) | Anywhere | — |
| `triangle` | `#` prefix | Top of each file (`.node`/`.ele`/`.poly`) | Yes |
| `ugrid` | None (fixed columnar numeric format) | n/a | — |
| `unv` | None general-purpose — dataset 2414's five ID-label records are a structural slot, not a comment | n/a | — |
| `vti` | XML `<!-- -->` | Anywhere in the document | Yes |
| `vtk` / `vtk42` / `vtk51` | The legacy format's title line (line 2, ≤256 chars) is the format's own free-text slot | Fixed line 2 | Yes |
| `vts` | XML `<!-- -->` | Anywhere in the document | Yes |
| `vtr` | XML `<!-- -->` | Anywhere in the document | Yes |
| `vtkhdf` | The `VTKHDF` group's `meshioplusplus:provenance` attribute (plain UTF-8 text) | Group metadata; ignored by `vtkHDFReader` | Yes |
| `vtm` | XML `<!-- -->` (the index file only; each piece carries its own) | Anywhere in the document | Yes |
| `pvd` / `pvtu` / `pvtp` | XML `<!-- -->` (the index file only; each piece carries its own) | Anywhere in the document | Yes |
| `radioss` | `#` or `$` prefix (`#include` and `#enddata` excepted) | Anywhere; the writer puts it after `#RADIOSS STARTER` (v16.17.0) | Yes |
| `vtp` | XML `<!-- -->` | Anywhere in the document | Yes |
| `vtu` | XML `<!-- -->` | Anywhere in the document | Yes |
| `wkt` | None (the OGC WKT grammar has no comment token) | n/a | — |
| `xdmf` | XML `<!-- -->` | Anywhere in the document | — |
| `xyz` | `#` prefix | Top of file, before the `# x y z ...` column header | Yes |
| `z88` | None (text after the header line's integers is ignored: Z88's own tools write a comment there) | The end of `z88i1.txt`'s header line | Yes |
| `zarr` | The root group's `meshioplusplus:provenance` attribute (plain JSON) | Store metadata; recovered without importing zarr | Yes |

Two formats carry a related but **structurally distinct** record that this table does not count as "the tag": `med`'s `DES` mesh-description field defaults to `"Mesh created with meshio++"` (both engines agree; user-overridable, so it is data, not a fixed credit) and `unv`'s dataset-2414 field-header records always read `meshioplusplus` on five fixed ID lines (a label field the format requires, not a comment).

`nastran` is the one documented exception to "both engines emit character-identical bytes": the C++ writer emits a literal sentinel comment (`meshioplusplus-cpp-nastran`) first and the provenance tag as a second `$` line after it; the Python writer emits only the tag. The sentinel is vestigial since v16.1.0 — the C++ reader of earlier releases accepted only files carrying it, and today's reads any deck — see [`nastran.md`](./formats/nastran.md).

---

## Native acceleration and fallbacks

meshio++ ships a C++ core (`meshioplusplus._core`, built with pybind11 + scikit-build-core). Most formats read and write through the C++ core with zero-copy numpy at the I/O boundary; each has a pure-Python fallback that is used automatically when the C++ path can't handle a file or when the extension was built without an optional dependency:

- **HDF5** (`cgns`, `h5m`, `hmf`, `med`, `vtkhdf`, and XDMF `data_format="HDF"`) — C++ when built with `MESHIOPLUSPLUS_WITH_HDF5`, otherwise `h5py`. For `med`, the C++ core covers the mesh-representation part (points, tags, families — including ones synthesized from named regions or `gmsh:physical`, same-type block consolidation — metadata, node orientation, `POG` ragged polygons) and defers the field/bitmask/multi-mesh constructs to the Python reference; see [`med.md`](./formats/med.md#quirks-limitations). `cgns` is a genuine CGNS/SIDS-compliant subset since v9.8.0 (readable by cgnslib/ParaView/VTK), covering the fixed-node-count element types and, since v9.21.0, polyhedral `NGON_n`/`NFACE_n` sections in both directions; see [`cgns.md`](./formats/cgns.md).
- **netCDF** (`exodus`) — C++ when built with `MESHIOPLUSPLUS_WITH_NETCDF`, otherwise `netCDF4`.
- **ADIOS2** (`vtx`) — C++ when built with `MESHIOPLUSPLUS_WITH_ADIOS2` (never in the release wheels), otherwise the `adios2` package; a core built with ADIOS2 does not fall back, since the package would load a second copy of the ADIOS2 libraries.
- **TecIO** (`szplt`) — C++ when built with `MESHIOPLUSPLUS_WITH_TECIO`, otherwise a shared TecIO through ctypes (`MESHIOPLUSPLUS_TECIO_LIBRARY`).
- **zlib** (VTU zlib compression) — C++ when built with `MESHIOPLUSPLUS_WITH_ZLIB`, otherwise the Python stdlib.

`mdpa` is the one format where the Python API deliberately does **not** prefer the C++ core for reading: only the pure-Python reference produces MDPA's `mesh.misc_data`, `mesh.geometries_block` and nested-by-cell-type `cell_data`. The C++ reader/writer exists (and is what the C API / Fortran / Julia / R / WebAssembly / native CLI use), and `mdpa.write` does use it for meshes carrying none of those extras — see [MDPA](./formats/mdpa.md#c-core).

Behaviour and file compatibility are identical either way; the native paths are only faster. Install the optional runtime deps with `pip install meshioplusplus[all]`.

### When the native path declines

A format's Python shim asks `meshioplusplus._fallback.core_declined` what to do with anything the C++ core raises, so a fallback is never silent and never unconditional:

| Raised by the core | Meaning | Behaviour |
| --- | --- | --- |
| `ReadError` / `WriteError` | A recognised decline: a malformed file, or a construct the core deliberately does not handle | Logged at `DEBUG`, falls back to the Python reference |
| `TypeError`, `MemoryError`, `RecursionError` | Never a decline: a stale build, a genuine user error, or retrying a huge file on the more memory-hungry twin | Propagates |
| Anything else (`RuntimeError`, `AttributeError`, …) | The fast path broke on something it did not classify: a bug in the core, since v16.14.0 | Logged at `WARNING` naming the format, path and cause, then falls back |

Since v16.14.0 every native reader entry point (the Python `_core.*_read` bindings and the registry every flat binding reads through) rethrows a parser's `std::` exception (`std::stoll` on a bad token, `.at()` past an end) as a `ReadError` naming the format (`detail/read_guard.hpp`), so a malformed file is always the first row and the last row is left for genuine defects. The fuzz campaign holds the readers to that contract ([fuzzing](./fuzzing.md)).

The fallback is also skipped where the Python twin cannot answer the same question: a non-default `time_step` (Gmsh, XDMF, MED, EnSight, CGNS) or an OpenFOAM `region` re-raises the core's error rather than quietly returning step 0 or a single region.

To see the declines, enable logging (`logging.basicConfig(level=logging.DEBUG)` shows the `meshioplusplus` logger). To prove a file was really handled by the core, set `MESHIOPLUSPLUS_STRICT_CORE=1` (also `on`, `true`, `yes`): every decline then re-raises instead of falling back, and is logged at `WARNING`. The ambiguous-extension loop in `read()` still moves on to the next candidate format, since that is about format ambiguity rather than core versus Python. Since v16.20.0 that loop tries first the candidate `sniff_format` recognises from the file's first 512 bytes (a Gmsh `.msh` goes to `gmsh` before `ansys` and `freefem`), still trying every candidate; and a reader that can tell from a file's first bytes that it is not its format -- a VTK XML file with a compressor the core cannot decode, a Gmsh file with `$Periodic`, a `.msh` that is not Fluent -- declines before parsing the rest. The operations (`clean`, `smooth`, …) have their own contract, `core_op_declined`, because an operation's Python twin is a different algorithm rather than the same file read differently:

| Raised by the core | Meaning | Behaviour |
| --- | --- | --- |
| `ImportError` (no compiled core) | A pure-Python install | The twin runs, silently |
| `NotImplementedError` (C++ `Unsupported`), or `ReadError`/`WriteError` from handing the mesh to the core | An input the native kernel deliberately leaves to the twin, or one the numpy conversion cannot carry (a ragged block, a string array) | Logged at `DEBUG`, the twin runs |
| `ValueError`, `TypeError`, `MemoryError`, `RecursionError` | A bad argument | Propagates: the twin never answers a question the core rejected |
| Anything else | A bug in the core | Logged at `WARNING`, the twin runs; re-raised under `MESHIOPLUSPLUS_STRICT_CORE=1` |

---

## Format-specific write options

All writers are called as `meshioplusplus.write(filename, mesh, file_format=..., **kwargs)` or `mesh.write(filename, **kwargs)`. The `**kwargs` depend on the format.

### Gmsh (`.msh`)

```python
meshioplusplus.gmsh.write(filename, mesh,
    fmt_version="4.1",   # "2.2", "4.0", or "4.1"
    binary=True,
    float_fmt=".16e",
)
```

Use `file_format="gmsh22"` to write version 2.2 via the generic `meshioplusplus.write`.

### VTU (`.vtu`)

```python
meshioplusplus.vtu.write(filename, mesh,
    binary=True,
    compression="zlib",   # "zlib", "lzma", or None
    header_type=None,     # "UInt32" or "UInt64"
)
```

### VTI (`.vti`)

```python
meshioplusplus.vti.write(filename, mesh,   # mesh must be a dense lattice
    binary=True,
    compression="zlib",   # "zlib", "lz4", "zstd", or None
    header_type=None,     # "UInt32" or "UInt64"
)
```

### VTS (`.vts`)

```python
meshioplusplus.vts.write(filename, mesh,   # mesh must be a dense lattice (write only -- read accepts any structured grid)
    binary=True,
    compression="zlib",   # "zlib", "lz4", "zstd", or None
    header_type=None,     # "UInt32" or "UInt64"
)
```

### VTR (`.vtr`)

```python
meshioplusplus.vtr.write(filename, mesh,   # mesh must be a UNIFORM dense lattice (write only -- read accepts any monotonic per-axis spacing)
    binary=True,
    compression="zlib",   # "zlib", "lz4", "zstd", or None
    header_type=None,     # "UInt32" or "UInt64"
)
```

### VTM (`.vtm`)

```python
meshioplusplus.vtm.write(filename, mesh,   # any mesh with one or more cell blocks -- no lattice restriction
    binary=True,
    compression="zlib",   # "zlib", "lz4", "zstd", or None
    header_type=None,     # "UInt32" or "UInt64"
)
```

### PVTU (`.pvtu`) and PVTP (`.pvtp`)

```python
meshioplusplus.pvtu.write(filename, mesh,   # carved by the integer cell_data array `part_key`, else one piece
    binary=True,
    compression="zlib",   # "zlib", "lz4", "zstd", or None
    header_type=None,     # "UInt32" or "UInt64"
    part_key="partition:part",
)
meshioplusplus.pvtu.write_pieces(filename, pieces)  # e.g. partition(mesh, n, ghost_layers=1); every piece must declare identical arrays

meshioplusplus.pvtu.read(filename, piece=None, ghosts="keep")  # piece=k keeps one; ghosts="drop" removes vtkGhostType cells (also meshioplusplus.read(..., ghosts="drop"))
```

`meshioplusplus.pvtp` takes the same arguments over `.vtp` pieces.

### PVD (`.pvd`)

```python
meshioplusplus.pvd.write(filename, mesh, binary=True, compression="zlib")  # a one-step collection
meshioplusplus.write_sequence("run.pvd", steps)                            # many steps, streamed

meshioplusplus.read(filename, time_step=-1, piece=None)  # a step (timestep=), and/or one part (part=) of it
```

### VTKHDF (`.vtkhdf`)

```python
meshioplusplus.vtkhdf.write(filename, mesh,
    compression="gzip",             # "gzip" or None; datasets under 4 KiB are never compressed
    compression_opts=4,             # gzip level 0-9
    dataset_type="UnstructuredGrid",  # "PolyData", "PartitionedDataSetCollection" or "MultiBlockDataSet"
    version=None,                   # None = oldest covering version, or (major, minor)
)

meshioplusplus.read(filename, time_step=-1, piece=None)  # a step, and/or one piece
```

### VTK (`.vtk`)

```python
meshioplusplus.vtk.write(filename, mesh,
    binary=True,
    # For version selection use file_format="vtk42" or "vtk51"
)
```

`file_format="vtk"` writes VTK 5.1. `file_format="vtk42"` or `"vtk51"` select specific versions.

### XDMF (`.xdmf`, `.xmf`)

```python
meshioplusplus.xdmf.write(filename, mesh,
    data_format="HDF",        # "HDF", "XML", or "Binary"
    compression="gzip",       # h5py compression filter (HDF only)
    compression_opts=4,       # compression level
)
```

With `data_format="HDF"`, meshio++ writes a companion `.h5` file alongside the `.xdmf`. With `"XML"`, all data is embedded in the XML. With `"Binary"`, data is written to separate `.bin` files.

### Medit (`.mesh`)

```python
meshioplusplus.medit.write(filename, mesh,
    float_fmt=".16e",
)
```

### PLY (`.ply`)

```python
meshioplusplus.ply.write(filename, mesh,
    binary=True,
)
```

### STL (`.stl`)

```python
meshioplusplus.stl.write(filename, mesh,
    binary=False,
)
```

### MED (`.med`)

```python
meshioplusplus.med.write(filename, mesh,
    med_version="4.1.0",   # MAJ.MIN.REL written to INFOS_GENERALES
)
```

MED does not support compression. `meshioplusplus.med.read_med_multi`/ `write_med_multi` read/write files containing several meshes — see [`med.md`](./formats/med.md). Since v9.6.0 MED is also a Phase-1 [named region](./regions.md) format (`FAS`/`GRO` group names ↔ `Point`/`Cell` regions, no side regions), carries the optional `NUM` global numbering as `point_data`/`cell_data["med:num"]`, and rejects a file written by a newer MED major version with a named error.

### AnsysInp (`.cdb`, `.inp`)

`meshioplusplus.ansysInp.read(filename, lenient=False)` / `meshioplusplus.ansysInp.write(filename, mesh)`: `lenient` skips elements whose type has no meshio++ cell instead of failing. See the [`.inp` note](#format-table) above for the Abaqus extension collision.

### OpenFOAM (`.foam`)

`meshioplusplus.openfoam.read(filename)` / `meshioplusplus.openfoam.write(filename, mesh, binary=False, label_bits=32, scalar_bits=64)`. `write` creates `<case>/constant/polyMesh/`; it is the only meshio++ writer that produces a directory, and needs the compiled core (there is no Python fallback writer). See [OpenFOAM](formats/openfoam.md#binary-write).

### glTF (`.glb`, `.gltf`)

```python
meshioplusplus.gltf.write(filename, mesh,
    split_angle=30.0,       # degrees in [0, 180]: crease angle above which points are duplicated
    normal_weight="angle",  # "angle" or "area"
    normals=True,           # False: no NORMAL, shared vertices
    fields=True,            # point_data (1-4 components) as raw _NAME attributes
    color_by=None,          # point_data / cell_data array baked into COLOR_0 (linear, unlit material)
    component=None, cmap="viridis", vmin=None, vmax=None, nan_color="#808080", unlit=True,
    up_axis="auto",         # "auto", "x", "y", "z": the source axis that points up
    recenter=True,          # subtract the bounding-box centre (carried in the root node)
    scale=1.0,              # source unit -> metres (root node scale)
    by_region=True,         # one node per cell region
    container="auto",       # "auto" (suffix), "glb" or "gltf"
)
```

### PCD (`.pcd`)

```python
meshioplusplus.pcd.write(filename, mesh,
    data=None,             # "ascii", "binary" (default) or "binary_compressed"
    point_dtype="float32", # "float32" (PCL's), "float64" or "keep"
    binary=True,           # picks binary/ascii when data is omitted
)
meshioplusplus.pcd.read(filename, drop_invalid=False)  # drop the points with a non-finite coordinate
```

### XYZ (`.xyz`, `.xyzn`, `.xyzrgb`, `.asc`, `.pts`, `.txt`)

```python
meshioplusplus.xyz.write(filename, mesh,
    float_fmt=None,   # e.g. ".16e"; default ".9g" for float32 columns, ".17g" otherwise
)
meshioplusplus.xyz.read(filename, columns=None, delimiter=None)  # columns: ["x", "y", "z", "nx", ...]
```

### CGNS (`.cgns`)

```python
meshioplusplus.cgns.write(filename, mesh,
    compression="gzip",
    compression_opts=4,
)
```

### Nastran (`.bdf`)

```python
meshioplusplus.nastran.write(filename, mesh,
    point_format="fixed-large",   # or "fixed-small", "free"
    cell_format="fixed-small",
)
```

### FLAC3D (`.f3grid`)

```python
meshioplusplus.flac3d.write(filename, mesh,
    float_fmt=".16e",
    binary=False,
)
```

`ZGROUP`/`FGROUP` cell groups round-trip as named [regions](./regions.md) called `<zone|face>:<name>:<slot>` — see [FLAC3D](./formats/flac3d.md#data-mapping) for the naming rule and for the `cell_sets` index convention, which changed in v10.36.0.

### SU2 (`.su2`)

`meshioplusplus.su2.write(filename, mesh)` — no extra options. A mesh whose `cell_data["su2:zone"]` carries more than one distinct value writes a single-file [multizone](./formats/su2.md#multizone-nzone) mesh (`NZONE=`/`IZONE=`, v15.5.0); a named boundary marker round-trips as a [region](./regions.md) instead of collapsing into a bare `su2:tag` id.

### AVS-UCD (`.avs`)

`meshioplusplus.avsucd.write(filename, mesh)` — no extra options.

### Elmer (mesh directory)

`meshioplusplus.elmer.write(dirname, mesh)` — no extra options; a directory has no extension, so `meshioplusplus.write` needs `file_format="elmer"`. Cell regions become bodies and boundaries, side regions boundary elements, with parents regenerated; `cell_data["partition:part"]` adds a `partitioning.N` directory; see [`elmer.md`](./formats/elmer.md#writing).

### FEBio (`.feb`)

`meshioplusplus.febio.write(filename, mesh)` — no extra options. Spec 4.0 with a placeholder material per domain; 2-D blocks on solid faces become `<Surface>`s, line blocks `<Edge>`s or `<DiscreteSet>`s, regions node sets, element sets and surfaces; point and cell data `<MeshData>`. See [`febio.md`](./formats/febio.md#writing).

### Code_Aster (`.mail`)

`meshioplusplus.code_aster.write(filename, mesh)` — no extra options. Point and cell regions become `GROUP_NO`/`GROUP_MA`, with names sanitised to 24 letters, digits and `_`; every line stays within the 80 columns Code_Aster reads; see [`code_aster.md`](./formats/code_aster.md#writing).

### Patran neutral (`.pat`, `.out`)

`meshioplusplus.patran.write(filename, mesh)` — no extra options. Point and cell regions become named components (one per name), the element property comes from `patran:property`, and coordinates keep ten significant digits (`E16.9`); see [`patran.md`](./formats/patran.md#writing).

### Femap neutral (`.neu`)

`meshioplusplus.femap.write(filename, mesh)` — no extra options. The Femap 8.2 layout: properties from `femap:property`/`femap:type` and their regions, other regions as groups; results are not written. See [`femap.md`](./formats/femap.md#writing).

### MFEM (`.mesh`)

```python
meshioplusplus.mfem.write(filename, mesh, grid_functions=False)
```

`.mesh` defaults to Medit, so pass `file_format="mfem"` to `meshioplusplus.write`. `grid_functions=True` also writes each data array as `<stem>.<name>.gf`. Quadratic cells become an `H1_<d>D_P2` nodes space and named cell regions v1.3 attribute sets; see [`mfem.md`](./formats/mfem.md#writing).

### Z88 (`z88i1.txt`)

`meshioplusplus.z88.write(filename, mesh, stubs=False)` — the Z88OS v15 structure file. The format follows the file name `z88i1.txt`; any other name needs `file_format="z88"`. Element types come from `z88:type`, else from the cell type; the `z88:` deck arrays (constraints, materials, element parameters, integration orders, surface loads) are written back beside it, and cell and point regions as Z88Aurora's `z88sets.txt`; other data arrays and side regions are dropped. `stubs=True` also writes empty `z88i2.txt` and `z88i5.txt` when there is nothing to put in them. See [`z88.md`](./formats/z88.md#writing).

### Abaqus (`.inp`)

Abaqus is one of the three Phase-1 [named region](./regions.md) formats (with gmsh and MED), and the only one that can express a **side set**: `*NSET` → point regions, `*ELSET` → cell regions and `*SURFACE, TYPE=ELEMENT` → side regions, in both the C++ core and the Python reference. Abaqus names its groups but has no integer id for them, so a region's `tag` is not preserved. Face identifiers (`S1`..`S6`) are remapped to meshio++'s own facet numbering — the two differ, and the per-type table is spelled out in [Named regions](./regions.md#abaqus-face-identifiers).

`meshioplusplus.abaqus.write(filename, mesh)` — no extra options.

### LS-DYNA (`.k`, `.key`, `.dyn`)

`meshioplusplus.lsdyna.write(filename, mesh)` — no extra options. Cell regions with a dimension become `*PART` cards (the tag is the `pid`), every other region becomes a `*SET_*_LIST`, and each part gets placeholder section and material ids; see [`lsdyna.md`](./formats/lsdyna.md#writing).

### DOLFIN-XML (`.xml`)

`meshioplusplus.dolfin.write(filename, mesh)` — no extra options. Both `point_data` (since v9.9.0, as `dim="0"` mesh functions) and `cell_data` are written to sibling `<stem>_<name>.xml` files; see [`dolfin.md`](./formats/dolfin.md#file-structure).

### GiD (`.post.msh` / `.post.res`, `.post.bin`, `.post.h5`)

```python
meshioplusplus.gid.write(filename, mesh,
    mode="auto",             # "auto", "ascii", "binary", or "hdf5"
    analysis_name="meshio++",
    step=1.0,
)
```

See [`gid.md`](./formats/gid.md). Reading takes only `time_step` (selecting one step of a multi-step results file).

---

## CLI format names

When using `meshioplusplus convert -o <format>`, use one of the format names from the first column of the table above (e.g. `gmsh`, `gmsh22`, `vtk`, `vtk42`, `vtu`, `xdmf`, …).
