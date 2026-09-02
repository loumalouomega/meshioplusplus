"""Headless Blender demo: read a mesh with meshio++ and inspect it in Blender.

    blender --background --python example/blender/import_demo.py -- example/example.msh

A script rather than a notebook, deliberately: the notebooks in ``example/`` are
re-executed in CI, and ``bpy`` is a several-hundred-megabyte dependency pinned
to one CPython minor. This is the ``example/physicsnemo/`` precedent — a worked
directory that is run on purpose, not on every docs build.
"""

import sys

import bpy

import meshioplusplus as mio


def main(path):
    mesh = mio.read(path)
    print(
        f"read {path}: {len(mesh.points)} points, "
        f"{sum(len(b) for b in mesh.cells)} cells in {len(mesh.cells)} block(s)"
    )

    data = mio.to_blender(mesh, name="demo")
    obj = bpy.data.objects.new(data.name, data)
    bpy.context.collection.objects.link(obj)

    print(
        f"blender: {len(data.vertices)} verts, {len(data.polygons)} polygons, "
        f"{len(data.edges)} edges"
    )
    sizes = sorted({p.loop_total for p in data.polygons})
    print(f"polygon sizes present: {sizes}  (n-gons are kept, not triangulated)")
    for layer in data.attributes:
        print(f"  attribute {layer.name!r}: {layer.data_type} on {layer.domain}")

    # Blender's own topology checker: True means it had to repair something.
    print("validate() had to fix something:", data.validate(verbose=False))

    back = mio.from_blender(data)
    print(f"back out: {[(b.type, len(b)) for b in back.cells]}")


if __name__ == "__main__":
    args = sys.argv[sys.argv.index("--") + 1 :] if "--" in sys.argv else []
    if not args:
        raise SystemExit("usage: blender --background --python import_demo.py -- MESH")
    main(args[0])
