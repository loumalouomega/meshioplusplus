# Skin to CAD

A triangle skin — an STL from a scan or a 3-D-printing download, the output of [`extract_skin`](/extract_skin) or [`isosurface`](/isosurface) — is thousands of flat facets nobody can select, dimension or edit in a CAD tool. Reverse engineering it means recovering the surfaces it was made of: planes, cylinders, cones, spheres and tori, then freeform patches, sewn into a solid. meshio++ does the mesh-side half of that work with no CAD kernel in the loop:

1. [**`segment`**](/segment) groups the facets into smooth patches bounded by the feature edges, with a curvature test that separates a fillet from its tangent planes;
2. [**`fit_primitives`**](/fit_primitives) fits a plane, sphere, cylinder, cone or torus to each patch and reports the residual, refusing any fit outside the tolerance;
3. the [**`step`**](/formats/step) writer hands the closed skin to CAD as a faceted solid.

```python
import meshioplusplus as mp

skin = mp.repair(mp.read("part.stl"))                 # closed and consistently wound
seg = mp.segment(skin, curvature_tolerance=1.0, min_patch_cells=20)
fit, report = mp.fit_primitives(seg, tolerance=2e-3, return_report=True)

fd = fit.field_data
for pid, t, p in zip(fd["fit:patch_id"], fd["fit:patch_type"], fd["fit:patch_params"]):
    if t == 3:                                          # a cylinder: a hole or a boss
        print(f"patch {pid}: radius {p[6]:.3f} along {p[3:6]}")

mp.write("part.step", skin)                            # opens in CAD as a solid
```

![A noisy turned part segmented into six patches: base, wall, fillet, shoulder, cone and dome](/images/skin_to_cad_patches.png)

The part above is a turned part tessellated into 2,240 facets with every point moved by up to `1e-3`. Its base, shoulder, cone and dome meet their neighbours at creases, but the toroidal fillet is tangent to the wall and to the shoulder: the feature edges alone find four patches, and the curvature test finds the six. Each patch then fits its primitive within `2e-3` of its size — the wall a cylinder of radius 2.0000, the fillet a torus of radii 1.5001 and 0.4998, the cone a half-angle of 21.80 degrees, the dome a sphere of radius 1.5006:

![The same part coloured by fitted type: two planes, a cylinder, a torus, a cone and a sphere](/images/skin_to_cad_fits.png)

The same chain is a [pipeline](/pipeline) — `[{"Op": "Segment"}, {"Op": "FitPrimitives"}]` — and three verbs of both CLIs (`segment`, `fit-primitives`, and `convert` to a `.step` file); every language surface and the MCP server expose it. Notebook [34](https://github.com/loumalouomega/meshioplusplus/blob/main/example/python/34_skin_to_cad.ipynb) builds that part, segments it, fits it and writes it as STEP.

## What is not done yet

Writing the fitted surfaces themselves as CAD faces — a cylinder face where the skin has a ring of facets — needs a CAD kernel to trim, sew and validate them, and so does reading STEP. Both wait on one decision, whether to add an optional OpenCASCADE backend; the [findings](/spikes#the-opencascade-dependency) record where it stands. Rolling-ball fillets, swept fillets and threads (a second pass over the fitted patches) and a B-spline fallback for what no primitive explains are on the [roadmap](/roadmap) behind it. Until then the solid's faces stay flat, and the fits are measurements.

The reference for the whole pipeline is [stl2curves](https://github.com/Duvdivan/stl2curves) (MIT), a CLI and FreeCAD add-on built on OpenCASCADE; meshio++ read its README only, and no code was taken from it.
