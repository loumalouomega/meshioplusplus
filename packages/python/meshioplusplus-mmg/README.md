# meshioplusplus-mmg

The MMG add-on for [meshio++](https://github.com/loumalouomega/meshioplusplus). It registers a remesher with the base wheel's shared core, which gives `adapt`, `discretize_levelset` and `mesh_domain` their MMG backend. `move` needs an MMG built with ELAS, which makes the combination GPL-3.0, so this wheel is built without it and `move` raises naming that.

```
pip install "meshioplusplus[mmg]"
```

The add-on is pinned to the exact base version and is licensed `MIT AND LGPL-3.0-or-later`; libmmg is linked dynamically and shipped as a shared library, so a recipient can replace it. See `doc/python_packages_plan.md` in the repository.
