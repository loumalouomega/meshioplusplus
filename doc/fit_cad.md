# Fit CAD

`fit_cad(mesh, path)` turns the fitted patches of a closed triangle skin into a CAD solid: every patch that [`fit_primitives`](/fit_primitives) accepted becomes a plane, sphere, cylinder, cone or torus face trimmed by the patch boundary, the faces are sewn into a solid that OpenCASCADE checks, and the solid is written as STEP. It is the last step of [skin to CAD](/cad_fit), after [`segment`](/segment) and `fit_primitives`, and it needs the optional [OpenCASCADE backend](#getting-opencascade).

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

- **Faces.** A patch whose fit was accepted, and whose type is in `types`, becomes one face on its fitted surface. Every other surface cell (an unfitted patch, a type left out, a cell with no patch) becomes planar facets, one per cell, as the faceted writer writes them.
- **Edges.** Where two patches meet, the mesh edges between them form chains from one junction (a point where three or more patches meet) to the next; a closed chain is split in two. Between two analytic faces the edge is the branch of the two surfaces' exact intersection that lies within the tolerance of the chain, when `exact_edges` is set and one does (a plane meeting a cylinder gives a circle). Otherwise, notably where a fillet meets its wall tangentially and no intersection is stable, the edge is a B-spline through the chain, smoothed and moved onto both surfaces, kept only if it stays within the tolerance of both. Beside a facet, the edge is the mesh's own straight edge.
- **Vertices.** A junction among analytic faces only is moved onto all of their surfaces; one beside a facet keeps the mesh position.
- **Solid.** Each face gets parameter-space curves that are continuous around its loops (so a boundary running along a cylinder's seam or a torus's equator does not fall apart), OpenCASCADE adds the seams and repairs the faces, sews them and closes each shell into a solid, one per edge-connected piece of the skin. A patch with no boundary at all, a whole sphere or torus, is its surface's natural face.

## Checks and demotion

A face is not kept when OpenCASCADE's B-rep check rejects it, when sewing leaves an edge of it free, when its boundary cannot be built, or when it strays: the solid is re-meshed and a face with a point further than `tolerance` from the skin, or a patch with a skin point further than that from the solid, fails. A failed patch is **demoted** to facets and the solid is rebuilt; after `max_demotion_rounds` rounds (default 8) every remaining patch is faceted, which is always a valid solid for a closed skin. So a feature the fitter cannot explain, or a face the kernel cannot trim, comes back faceted, never silently wrong, and `num_demoted` counts the features dropped. `MESHIOPLUSPLUS_LOG_LEVEL=debug` logs why each patch was demoted.

`tolerance` defaults to three times the largest `fit:patch_tolerance` of an accepted patch, so it follows the scale `fit_primitives` worked at; with no accepted patch it is `1e-3` of the bounding-box diagonal. The re-meshing uses a chordal deflection of a quarter of the tolerance unless `mesh_deflection` says otherwise.

## Output and report

The return value is the solid re-meshed by OpenCASCADE as one `triangle` block, with points welded where faces share them exactly, and Int64 cell data `fit_cad:face` (the B-rep face), `fit_cad:type` (0 a facet, else the `fit:type` value: 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus) and `fit_cad:patch` (the `segment:patch` label, -1 for a cell with no patch). With `return_report=True` the report holds:

| Key | Meaning |
| --- | --- |
| `valid` | the solid passes OpenCASCADE's B-rep check |
| `tolerance` | the tolerance applied |
| `volume`, `skin_volume`, `volume_ratio` | the solid's volume, the skin's enclosed volume and their ratio |
| `max_deviation` | the symmetric Hausdorff distance between the skin and the re-meshed solid |
| `num_faces`, `num_plane` … `num_torus`, `num_faceted_faces` | B-rep faces in all, by type, and planar facets |
| `num_patches`, `num_demoted` | patches in the fit table, and accepted ones demoted to facets |
| `num_exact_edges`, `num_approx_edges` | edges between analytic faces built as exact intersections and as B-splines |
| `num_solids` | solids written, one per piece of the skin |

The STEP file is AP214 in millimetres. Its header timestamp is pinned and its product is named `meshio++ fit_cad`, so the same input writes the same bytes; the bytes do depend on the OpenCASCADE version. The work is serial, so nothing depends on the thread count.

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
    [--mesh-deflection D] [--max-demotion-rounds N] [--mesh-out FILE] [--json] [--quiet]
```

A skin without `segment:patch` is segmented with `--curvature-tolerance` and `--min-patch-cells`, and fitted with `--fit-tolerance` (relative, as in `fit-primitives`), before the solid is built. `--mesh-out` also writes the re-meshed solid. Both CLIs exit with 2 when the solid fails the B-rep check, as `check` does for a failed gate.

```sh
meshioplusplus to-cad part.stl part.step --curvature-tolerance 1 --min-patch-cells 20 --fit-tolerance 2e-3
```

## Pipeline

```json
{"Op": "FitCad", "Path": "part.step", "Tolerance": 0.0, "ExactEdges": true,
 "Types": ["plane", "cylinder"], "MeshDeflection": 0.0, "MaxDemotionRounds": 8}
```

The step writes the solid to `Path` when given, and the mesh continues as the re-meshed solid; it reports `Valid`, `VolumeRatio`, `MaxDeviation`, `NumFaces`, `NumFacetedFaces` and `NumDemoted`. A `Tolerance` or `MeshDeflection` of 0 is automatic; an empty `Types` facets everything.

## Other languages

```c
mio_fit_cad_opts opts;
mio_fit_cad_opts_init(&opts);
opts.types = (1u << MIO_FIT_PLANE) | (1u << MIO_FIT_CYLINDER);
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

The MCP server's `fit_cad` tool takes the same options plus `curvature_tolerance`, `min_patch_cells` and `fit_tolerance`, and reports `occt_available`; without the backend it answers with the named error and `occt_available: false`. There is no WASM binding.

## Not supported

- **Fillets, blends and threads** as their own features: a constant-radius fillet that `segment` separates and `fit_primitives` fits as a torus becomes a toroidal face, but a swept or variable fillet and a thread come back faceted ([roadmap](/roadmap) 6.3.1).
- **Freeform patches**: what no primitive explains stays faceted; a B-spline surface fit is roadmap 6.3.2.
- **Reading STEP or IGES** is not done yet (roadmap 6.2.1), although the backend it needs is now in place.
