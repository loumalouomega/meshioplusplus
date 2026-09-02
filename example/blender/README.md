# Blender example

A headless run of the bridge, for checking an install end to end.

```sh
blender --background --python example/blender/import_demo.py -- ../example.msh
```

It reads the bracket with meshio++, converts it to a `bpy.types.Mesh`, prints
what arrived (vertex/polygon/edge counts, the polygon sizes present, every
attribute with its type and domain), runs Blender's own `validate()` as a
topology check, and converts back.

The interesting line is the polygon sizes: a mesh with quads shows a `4` there.
Unlike the trimesh bridge, nothing is triangulated — Blender holds n-gons.

Needs meshio++ installed into Blender's Python, or the extension zip installed;
see [`doc/blender.md`](../../doc/blender.md).
