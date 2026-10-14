# Skin to CAD

A triangle skin — an STL from a scan or a 3-D-printing download, the output of [`extract_skin`](/extract_skin) or [`isosurface`](/isosurface) — is thousands of flat facets nobody can select, dimension or edit in a CAD tool. Reverse engineering it means recovering the surfaces it was made of: planes, cylinders, cones, spheres and tori, then freeform patches, sewn into a solid. meshio++ does the mesh-side half of that work with no CAD kernel in the loop, and the CAD half through an optional OpenCASCADE backend:

1. [**`segment`**](/segment) groups the facets into smooth patches bounded by the feature edges, with a curvature test that separates a fillet from its tangent planes;
2. [**`fit_primitives`**](/fit_primitives) fits a plane, sphere, cylinder, cone or torus to each patch and reports the residual, refusing any fit outside the tolerance;
3. the [**`step`**](/formats/step) writer hands the closed skin to CAD as a faceted solid;
4. [**`fit_cad`**](/fit_cad), with OpenCASCADE, turns the fitted patches into analytic CAD faces trimmed by the patch boundaries, and a patch no primitive explains into a [freeform B-spline face](/fit_cad#freeform-faces), sews them into a checked solid and writes STEP a CAD tool can edit, or [Kratos CAD JSON](/formats/cad_json) for isogeometric analysis; a patch it cannot build comes back faceted.

```python
import meshioplusplus as mp

skin = mp.repair(mp.read("part.stl"))                 # closed and consistently wound
seg = mp.segment(skin, curvature_tolerance=1.0, min_patch_cells=20)
fit, report = mp.fit_primitives(seg, tolerance=2e-3, return_report=True)

fd = fit.field_data
for pid, t, p in zip(fd["fit:patch_id"], fd["fit:patch_type"], fd["fit:patch_params"]):
    if t == 3:                                          # a cylinder: a hole or a boss
        print(f"patch {pid}: radius {p[6]:.3f} along {p[3:6]}")

mp.write("part.step", skin)                            # opens in CAD as a faceted solid
solid, rep = mp.fit_cad(fit, "part_cad.step", return_report=True)  # analytic faces (OpenCASCADE)
```

![A noisy turned part segmented into six patches: base, wall, fillet, shoulder, cone and dome](/images/skin_to_cad_patches.png)

The part above is a turned part tessellated into 2,240 facets with every point moved by up to `1e-3`. Its base, shoulder, cone and dome meet their neighbours at creases, but the toroidal fillet is tangent to the wall and to the shoulder: the feature edges alone find four patches, and the curvature test finds the six. Each patch then fits its primitive within `2e-3` of its size — the wall a cylinder of radius 2.0000, the fillet a torus of radii 1.5001 and 0.4998, the cone a half-angle of 21.80 degrees, the dome a sphere of radius 1.5006:

![The same part coloured by fitted type: two planes, a cylinder, a torus, a cone and a sphere](/images/skin_to_cad_fits.png)

The same chain is a [pipeline](/pipeline) — `[{"Op": "Segment"}, {"Op": "FitPrimitives"}]` — and three verbs of both CLIs (`segment`, `fit-primitives`, and `convert` to a `.step` file); every language surface and the MCP server expose it. Notebook [34](https://github.com/loumalouomega/meshioplusplus/blob/main/example/python/34_skin_to_cad.ipynb) builds that part, segments it, fits it, writes it as STEP and, where the core carries OpenCASCADE, builds its analytic solid.

## Analytic faces

With a core built `-DMESHIOPLUSPLUS_WITH_OCCT=ON`, [`fit_cad`](/fit_cad) writes the fitted surfaces themselves: on the turned part above, the base and the shoulder become planes, the wall a cylinder, the fillet a torus, the cone a cone and the dome a sphere, six faces in one valid solid within `1e-2` of the skin. Where a fillet meets its wall tangentially the edge between them is a B-spline; elsewhere it is the exact circle where the two surfaces meet. The same chain runs as the pipeline step `FitCad` and as the `to-cad` verb of both CLIs, and every surface but WASM exposes it.

## Fillets, threads and curled patches

A patch no primitive explains is not given up on: `fit_cad` tries it as a rolling-ball fillet between its two analytic neighbours (a pipe around the spine the two surfaces meet along once each is offset by the radius, so a fillet swept along a curved edge works as well as a straight one, and, when no one radius explains it, a ball whose radius varies along the edge), as a helical thread face (a straight flank or any rounded profile, about the axis of a neighbouring cylinder or cone or of the screw the patch fits on its own), and otherwise as a freeform B-spline face; a patch no plane or unrolled cylinder parametrises is split into pieces that one does first. A patch between two fitted faces, one of them freeform, is then tried as a rolling-ball blend against them. See [fillets, blends and threads](/fit_cad#fillets-blends-and-threads). The same backend [reads STEP](/formats/step#reading) and [IGES](/formats/iges) models as their triangulation and NURBS B-rep. The [findings](/spikes#the-opencascade-dependency) record how the backend was decided.

The reference for the whole pipeline is [stl2curves](https://github.com/Duvdivan/stl2curves) (MIT), a CLI and FreeCAD add-on built on OpenCASCADE; meshio++ read its README only, and no code was taken from it.
