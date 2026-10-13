# Fit CAD

`fit_cad(mesh, path)` turns the fitted patches of a closed triangle skin into a CAD solid: every patch that [`fit_primitives`](/fit_primitives) accepted becomes a plane, sphere, cylinder, cone or torus face trimmed by the patch boundary, a patch no primitive explains a [freeform B-spline face](#freeform-faces), the faces are sewn into a solid that OpenCASCADE checks, and the solid is written as STEP. It is the last step of [skin to CAD](/cad_fit), after [`segment`](/segment) and `fit_primitives`, and it needs the optional [OpenCASCADE backend](#getting-opencascade).

```python
import meshioplusplus as mp

skin = mp.repair(mp.read("part.stl"))
seg = mp.segment(skin, curvature_tolerance=1.0, min_patch_cells=20)
fit = mp.fit_primitives(seg, tolerance=2e-3)

solid, report = mp.fit_cad(fit, "part.step", return_report=True)
report["valid"], report["num_cylinder"], report["num_demoted"]
report["volume_ratio"], report["max_deviation"]   # against the skin
```

![The turned part of skin to CAD as a CAD solid: two planes, a cylinder, a torus fillet, a cone and a spherical dome, re-meshed and coloured by face type](/images/fit_cad_solid.png)

Unlike the faceted [`step`](/formats/step) writer, whose solid is one flat face per triangle, the STEP file `fit_cad` writes carries `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`, `SPHERICAL_SURFACE` and `TOROIDAL_SURFACE` faces a CAD tool can select, dimension and edit. A skin without `segment:patch` or the `fit:patch_*` table is segmented and fitted first with the default options; a noisy scan wants a curvature test and a `min_patch_cells` in the tens, so run those two yourself, or pass them to the [CLI](#cli).

## How the solid is built

- **Faces.** A patch whose fit was accepted, and whose type is in `types`, becomes one face on its fitted surface; with `freeform` (the default), a patch no primitive explains becomes a [B-spline face](#freeform-faces). Every other surface cell (a type left out, a freeform fit that failed, a cell with no patch) becomes planar facets, one per cell, as the faceted writer writes them.
- **Edges.** Where two patches meet, the mesh edges between them form chains from one junction (a point where three or more patches meet) to the next; a closed chain is split in two. Between two analytic faces the edge is the branch of the two surfaces' exact intersection that lies within the tolerance of the chain, when `exact_edges` is set and one does (a plane meeting a cylinder gives a circle). Otherwise, notably where a fillet meets its wall tangentially and no intersection is stable, the edge is a B-spline through the chain, smoothed and moved onto both surfaces, kept only if it stays within the tolerance of both. Beside a facet, the edge is the mesh's own straight edge.
- **Vertices.** A junction among analytic faces only is moved onto all of their surfaces; one beside a facet keeps the mesh position.
- **Solid.** Each face gets parameter-space curves that are continuous around its loops (so a boundary running along a cylinder's seam or a torus's equator does not fall apart), OpenCASCADE adds the seams and repairs the faces, sews them and closes each shell into a solid, one per edge-connected piece of the skin. A patch with no boundary at all, a whole sphere or torus, is its surface's natural face.

## Freeform faces

A patch `fit_primitives` left unexplained (no accepted primitive, or one with non-finite parameters) gets a B-spline surface fitted to it, so a smooth blend, a dent or a styled panel becomes one CAD face instead of hundreds of facets:

1. **A parametrisation.** The patch is flattened over its best-fit plane (principal components) or, when that folds, unrolled around a cylinder: the axis is the direction its normals are most perpendicular to (the smallest eigenvector of Σ nnᵀ), the radius a least-squares circle fit of the points around it, the angle cut at the widest gap. A flattening whose triangles do not all keep one orientation is not a parametrisation and is not used; a patch that folds over both (a doubly curled one) stays faceted.
2. **A height field.** The offset from the plane or the cylinder is resampled on a grid over the flattened bounding box, interpolated inside the patch's triangles and continued from the nearest boundary point outside them, so the surface extends past the patch's trim.
3. **The fit.** `GeomAPI_PointsToBSplineSurface` fits the grid (iso-parametric, C2, degrees up to `freeform_degree`, default 3, valid 3–8). The grid starts near √n points a side and doubles up to `freeform_max_control` (default 64, valid 4–1024) until every patch vertex lies within half the tolerance and every facet centre within the tolerance; a patch that never gets there, or a skin too coarse for any smooth surface to pass within it, stays faceted.
4. **The face.** The surface is oriented by the facet normals and goes through the same trimming, sewing, check and [demotion](#checks-and-demotion) as an analytic face; its `fit_cad:type` is 6 (`kFitCadBsplineType`), and the report's `num_bspline` counts it.

`freeform=False` keeps every unexplained patch faceted.

## Checks and demotion

A face is not kept when OpenCASCADE's B-rep check rejects it, when sewing leaves an edge of it free, when its boundary cannot be built, or when it strays: the solid is re-meshed and a face with a point further than `tolerance` from the skin, or a patch with a skin point further than that from the solid, fails. A failed patch is **demoted** to facets and the solid is rebuilt; after `max_demotion_rounds` rounds (default 8) every remaining patch is faceted, which is always a valid solid for a closed skin. So a feature the fitter cannot explain, or a face the kernel cannot trim, comes back faceted, never silently wrong, and `num_demoted` counts the features dropped. `MESHIOPLUSPLUS_LOG_LEVEL=debug` logs why each patch was demoted.

`tolerance` defaults to three times the largest `fit:patch_tolerance` of an accepted patch, so it follows the scale `fit_primitives` worked at; with no accepted patch it is `1e-3` of the bounding-box diagonal. The re-meshing uses a chordal deflection of a quarter of the tolerance unless `mesh_deflection` says otherwise.

## Output and report

The return value is the solid re-meshed by OpenCASCADE as one `triangle` block, with points welded where faces share them exactly, and Int64 cell data `fit_cad:face` (the B-rep face), `fit_cad:type` (0 a facet, else the `fit:type` value: 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus; 6 a freeform B-spline) and `fit_cad:patch` (the `segment:patch` label, -1 for a cell with no patch). With `return_report=True` the report holds:

| Key | Meaning |
| --- | --- |
| `valid` | the solid passes OpenCASCADE's B-rep check |
| `tolerance` | the tolerance applied |
| `volume`, `skin_volume`, `volume_ratio` | the solid's volume, the skin's enclosed volume and their ratio |
| `max_deviation` | the symmetric Hausdorff distance between the skin and the re-meshed solid |
| `num_faces`, `num_plane` … `num_torus`, `num_bspline`, `num_faceted_faces` | B-rep faces in all, by type (freeform B-splines included), and planar facets |
| `num_patches`, `num_demoted` | patches in the fit table, and accepted ones demoted to facets |
| `num_exact_edges`, `num_approx_edges` | edges between analytic faces built as exact intersections and as B-splines |
| `num_solids` | solids written, one per piece of the skin |

The STEP file is AP214 in millimetres. Its header timestamp is pinned and its product is named `meshio++ fit_cad`, so the same input writes the same bytes; the bytes do depend on the OpenCASCADE version. The work is serial, so nothing depends on the thread count.

## The solid as NURBS

`solid.info["brep"]` (C++ `FitCadResult::mBrep`, a `BrepInfo`) holds the solid with every face converted to NURBS (`BRepBuilderAPI_NurbsConvert`): the surfaces, their trimming loops and the edges with the faces they bound, as a [STEP read](/formats/step#reading) gives them. Writing it as [Kratos CAD JSON](/formats/cad_json) (`mp.write("part.cad.json", solid)`, `to-cad --cad-json`) makes a scanned part Kratos IGA geometry.

## Getting OpenCASCADE

[Open CASCADE Technology](https://github.com/Open-Cascade-SAS/OCCT) is licensed under the GNU LGPL 2.1 with the additional exception in its `OCCT_LGPL_EXCEPTION.txt`. meshio++ stays MIT: the backend is **off by default and bring-your-own** (never vendored or downloaded), linked dynamically, and absent from the PyPI wheels, the WASM build and the release CLI binaries. A build that enables it links an LGPL library; distributing that build carries the LGPL's obligations for OpenCASCADE, not for meshio++. OpenCASCADE 7.6 or later is supported (7.6 is Ubuntu 24.04's package; the CI leg tests it).

```bash
sudo apt-get install libocct-foundation-dev libocct-modeling-data-dev \
    libocct-modeling-algorithms-dev libocct-data-exchange-dev libtbb-dev
CMAKE_ARGS="-DMESHIOPLUSPLUS_WITH_OCCT=ON" pip install meshioplusplus --no-binary meshioplusplus
```

From source, `build/configure.sh --with-occt` (or `-DMESHIOPLUSPLUS_WITH_OCCT=ON`) finds OpenCASCADE through its CMake config; point `OpenCASCADE_DIR` or `CMAKE_PREFIX_PATH` at an install outside the default prefixes, and make its libraries findable at run time as for any shared library. Debian's packages export toolkit targets that link TBB, hence `libtbb-dev`.

Without the backend, `fit_cad` raises `NotImplementedError` (C++ `Unsupported`) naming `-DMESHIOPLUSPLUS_WITH_OCCT=ON`; there is no pure-Python fallback, because trimming, sewing and checking surfaces is a CAD kernel's work. `meshioplusplus._fit_cad.has_occt()`, `fit_cad_has_occt()` in C++, `mio_has_occt()` in C and the language surfaces below say in advance whether a build has it.

## CLI

```sh
meshioplusplus to-cad IN OUT.step [--curvature-tolerance C] [--min-patch-cells N] \
    [--fit-tolerance T] [--tolerance T] [--no-exact-edges] [--types plane,cylinder,...] \
    [--mesh-deflection D] [--max-demotion-rounds N] [--no-freeform] [--freeform-degree P] \
    [--freeform-max-control N] [--mesh-out FILE] [--cad-json FILE] [--json] [--quiet]
```

A skin without `segment:patch` is segmented with `--curvature-tolerance` and `--min-patch-cells`, and fitted with `--fit-tolerance` (relative, as in `fit-primitives`), before the solid is built. `--mesh-out` also writes the re-meshed solid, `--cad-json` the solid as Kratos CAD JSON; `--no-freeform`, `--freeform-degree` and `--freeform-max-control` are the [freeform](#freeform-faces) options. Both CLIs exit with 2 when the solid fails the B-rep check, as `check` does for a failed gate.

```sh
meshioplusplus to-cad part.stl part.step --curvature-tolerance 1 --min-patch-cells 20 --fit-tolerance 2e-3
```

## Pipeline

```json
{"Op": "FitCad", "Path": "part.step", "Tolerance": 0.0, "ExactEdges": true,
 "Types": ["plane", "cylinder"], "MeshDeflection": 0.0, "MaxDemotionRounds": 8,
 "Freeform": true, "FreeformDegree": 3, "FreeformMaxControl": 64}
```

The step writes the solid to `Path` when given, and the mesh continues as the re-meshed solid; it reports `Valid`, `VolumeRatio`, `MaxDeviation`, `NumFaces`, `NumBspline`, `NumFacetedFaces` and `NumDemoted`. A `Tolerance` or `MeshDeflection` of 0 is automatic; an empty `Types` facets everything.

## Other languages

```c
mio_fit_cad_opts opts;
mio_fit_cad_opts_init(&opts);
opts.types = (1u << MIO_FIT_PLANE) | (1u << MIO_FIT_CYLINDER);
opts.freeform = 1;             /* the default; freeform_degree, freeform_max_control */
mio_fit_cad_report report;
if (mio_has_occt()) {
    mio_mesh* solid = mio_fit_cad(fitted, "part.step", &opts, &report);
}
```

```fortran
solid = fit%fit_cad(step_path='part.step', valid=ok, num_demoted=ndropped, stat=ierr)
```

```julia
has_occt() && (r = fit_cad(fit; step_path="part.step"); r.valid)
```

```r
if (mio_has_occt()) cad <- mio_fit_cad(fit, step_path = "part.step")
```

The C report carries `num_bspline`, as do Fortran's, Julia's and R's results (`freeform`, `freeform_degree` and `freeform_max_control` are their options). The MCP server's `fit_cad` tool takes the same options plus `curvature_tolerance`, `min_patch_cells`, `fit_tolerance` and `cad_json_output_path`, and reports `occt_available`; without the backend it answers with the named error and `occt_available: false`. There is no WASM binding.

## Not supported

- **Fillets, blends and threads** as their own features: a constant-radius fillet that `segment` separates and `fit_primitives` fits as a torus becomes a toroidal face, but a swept or variable fillet and a thread come back faceted ([roadmap](/roadmap) 6.3.1).
- **Doubly curled freeform patches**: a patch that folds over both its plane and its cylinder stays faceted rather than being split into flattenable pieces.
- **Fitting a STEP or IGES model**: those are [read](/formats/step#reading) as their own NURBS B-rep, not refitted.
