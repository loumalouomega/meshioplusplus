# Baking a field into a texture

`bake_texture` rasterizes a point or cell array into an RGBA image through the mesh's texture coordinates (`point_data["uv"]`, which [parameterize](./parameterize.md) makes), so a field travels as a picture: the glTF writer embeds it as the material's base-colour texture and the OBJ writer puts it beside the file with a `.mtl` (roadmap 2.4.8).

```python
import meshioplusplus

surface = meshioplusplus.parameterize(meshioplusplus.read("part.stl"))
surface.point_data["height"] = surface.points[:, 2]
textured, image = meshioplusplus.bake_texture(surface, "height", size=(1024, 1024), cmap="turbo")
meshioplusplus.write("part.glb", textured)   # the texture is embedded
meshioplusplus.write("part.obj", textured)   # part.obj, part.mtl, part.png
```

```sh
meshioplusplus bake-texture part_uv.vtu part.glb --array height --width 1024 --height 1024 --cmap turbo --image part.png
```

## What is drawn

- Each 2-D cell is drawn in texture space, `x = u * width` and `y = (1 - v) * height` (the image's top row is `v = 1`, as glTF and OBJ viewers sample it): a triangle as it is, a quad split along its 0-2 diagonal, a polygon fanned from its first node. Line and volume cells are not drawn.
- A pixel whose centre lies inside a triangle (all three barycentric weights non-negative) takes the array's value there: the barycentric blend of a point array, the cell's own value of a cell array. A later cell draws over an earlier one.
- A multi-component array reduces to `component`, else to its magnitude.
- The value maps through the colormap (the same tables as `colorize` and the glTF writer's `color_by`) over `vmin`..`vmax`, the finite range of the array by default. A pixel no cell covers, or whose value is not finite, stays transparent.

## Result

`(mesh, image)`: the input mesh with `field_data["texture:png"]` (the PNG bytes the glTF and OBJ writers embed) and the `(height, width, 4)` uint8 image; `return_report=True` adds `vmin`, `vmax`, `covered_pixels` and `warnings`. The C++ core (`operations/bake_texture.hpp`) and the numpy twin draw the same pixels.

## Bindings

Python, both CLIs (`bake-texture`) and the MCP server (`bake_texture`). The C API and the flat bindings (Fortran, Julia, R, Octave/MATLAB, WASM) and a pipeline step are not there yet.
