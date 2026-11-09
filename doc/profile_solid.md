# Profile solids

`profile_solid(sketch, features)` turns a closed [`fit_curves`](./fit_curves.md) sketch into a CAD solid through a CAD provider (roadmap 7.7.4): the sketch extruded along its plane's normal or revolved about an axis, then filleted, chamfered or shelled, written as STEP and tessellated. A turned or prismatic part becomes a clean solid from a few sections, and [`deviation`](./deviation.md) compares it to the scan.

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

section = meshioplusplus.slice(scan, origin, normal, chain=True)  # through the shaft's axis
sketch = meshioplusplus.fit_curves(section)
shaft, report = meshioplusplus.profile_solid(
    sketch, "revolve:0,0,0:0,0,1", "shaft.step", return_report=True
)
plate = meshioplusplus.profile_solid(outline, "extrude:0;fillet:0.5", scan=scan)  # to the scan's extent
meshioplusplus.deviation(scan, shaft)
```

## Design

The roadmap asked for a design pass before code: the smallest scripted API, and where it stops.

- **What it takes.** One sketch (a line mesh with `fit_curves`' `sketch:*` field data, from a section, a silhouette or a DXF file) and a list of features applied in order. The first makes the solid, `extrude` or `revolve`; the rest modify it: `fillet` and `chamfer` (every edge) and `shell` (hollow to a wall, removing the face the extrusion ended at). That covers the turned shaft, the prismatic plate and the cast cup the roadmap names, and nothing else.
- **Where it stops.** No sketch constraints, no feature history with rebuilds and rollback, no edge selection beyond "every edge", no patterns or threads: the interactive parametric modeller is a non-goal ([roadmap](./roadmap.md#non-goals-and-decisions-taken)). A part needing more goes to a CAD tool as STEP, which is what `profile_solid` writes.
- **Who builds it.** The core chains the sketch's pieces into closed loops in world coordinates (the outer loop the one of largest area, wound counter-clockwise about the sketch's normal, the holes clockwise), resolves an extrusion of 0 to the scan's extent along the normal, and hands that to a named profile-solid provider (`add_profile_solid_provider`, the `add_tet_mesher` pattern; `CadKernelProvider` is untouched). OpenCASCADE's provider, `"occt"`, builds the face, `BRepPrimAPI_MakePrism` or `BRepPrimAPI_MakeRevol`, `BRepFilletAPI_MakeFillet`, `BRepFilletAPI_MakeChamfer` and `BRepOffsetAPI_MakeThickSolid`, checks the solid, writes STEP (AP214, millimetres, a header free of time stamps, as [`fit_cad`](./fit_cad.md) writes) and tessellates it. The pinned OpenCASCADE gains the `TKFillet` and `TKOffset` toolkits for it. Without a provider the call throws `Unsupported` (`NotImplementedError` in Python) naming OpenCASCADE.
- **The features as text.** The CLIs, the pipeline and the language bindings take the features as one text, `;`-separated: `extrude:D` (`extrude:D:symmetric` half each way), `revolve:OX,OY,OZ:DX,DY,DZ[:ANGLE]` (360 by default), `fillet:R`, `chamfer:D`, `shell:W`. Python and C++ also take a list (`{"kind": "extrude", "distance": 10}`, `ProfileFeature`).

## The report

`volume`, `area`, `num_faces`, `valid` (OpenCASCADE's check of the solid), `provider` and `warnings`. A 2 × 1 rectangle extruded by 3 is a valid six-faced box of volume 6; filleted by 0.1 it loses volume and gains faces; shelled by 0.1 it is a cup of volume `6 - 1.8 · 0.8 · 2.9`; the half-section of a stepped shaft (radius 1 for 2, radius 0.6 for 1) revolved about its axis has volume `2.36 π`, and a scan of the shaft lies within 0.05 of it on more than 95% of its points (`tests/cpp/test_profile_solid.cpp`).

## Surfaces

C++ `profile_solid`, `profile_features_from_text`, `profile_feature_from_name` and the provider registry (`operations/profile_solid.hpp`); Python `meshioplusplus.profile_solid`, `profile_features_from_text` and `profile_solid_providers` (C++ core only); both CLIs (`profile-solid SKETCH OUT --features TEXT [--step PATH] [--scan SCAN] [--provider NAME] [--deflection F] [--json]`); the `ProfileSolid` [pipeline](./pipeline.md) step (`Features`, `Step`, `Scan`, `Provider`, `Deflection`); C `mio_profile_solid` (`mio_profile_feature`, 96 bytes, `mio_profile_feature_init`, `mio_profile_features_from_text`; `mio_profile_solid_report`, 64; `mio_profile_solid_warnings`, `mio_has_profile_solid_provider`); Fortran, Julia, R, Octave/MATLAB (features as text), WASM `profileSolid` (the full artifact, which carries OpenCASCADE); the MCP tool `profile_solid`.
