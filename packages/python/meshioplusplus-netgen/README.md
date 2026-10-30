# meshioplusplus-netgen

The Netgen add-on for [meshio++](https://github.com/loumalouomega/meshioplusplus). It registers a tetrahedral mesher with the base wheel's shared core, which gives `tetrahedralize` its Netgen backend: the surface's triangles are kept and the interior is filled. Without it, `tetrahedralize` fits the boundary with its native octree backend.

```
pip install "meshioplusplus[netgen]"
```

The add-on is pinned to the exact base version and is licensed `MIT AND LGPL-2.1-only`; Netgen's nglib is linked dynamically and shipped as shared libraries, so a recipient can replace it. The base wheel alone carries no copyleft code. See `doc/python_packages_plan.md` in the repository.
