# Catalyst 2 and Ascent

ParaView Catalyst 2 and Ascent hand an in situ pipeline a Conduit node holding a [Mesh Blueprint](../interop.md#conduit-mesh-blueprint). `contrib/catalyst/blueprint_series.py` reads it with `meshioplusplus.from_conduit` and writes one step per call, so a run lands as a file meshio++, ParaView and every registry writer read afterwards (roadmap 10.2). It needs `meshioplusplus` (and `h5py` for `.vtkhdf`) in the Python the pipeline runs in; Conduit comes with the host.

## Catalyst 2

The Catalyst API calls `catalyst_initialize`, `catalyst_execute` and `catalyst_finalize` with a node in the Catalyst protocol: the time in `catalyst/state/time`, the step in `catalyst/state/timestep`, and each channel's Blueprint in `catalyst/channels/<name>/data` (`type` `mesh` or `multimesh`). The script's three functions of those names take that node, from a Catalyst implementation or a simulation's Python adaptor:

```python
import blueprint_series as mio_catalyst

mio_catalyst.catalyst_initialize(node)   # node["catalyst/meshioplusplus/output"] = "run.vtkhdf"
for step in range(n):
    mio_catalyst.catalyst_execute(node)  # catalyst/state/time, catalyst/channels/grid/data
mio_catalyst.catalyst_finalize(node)
```

The output path is `catalyst/meshioplusplus/output` in the initialize node, else `$MIO_CATALYST_OUTPUT`, else `catalyst.vtkhdf`; `catalyst/meshioplusplus/channel` picks the channel (the first mesh channel by default). A `.vtkhdf` path is one transient [VTKHDF](../formats/vtkhdf.md) file, which holds one grid: a mesh that changes between steps is refused by name. A path with `{step}` writes one file per step, in any format meshio++ writes, numbered by `catalyst/state/timestep`.

## Ascent

As an Ascent Python extract (`type: python` with `params/file` naming the script), `ascent_data()` is the published Blueprint and `ascent_extract()` writes it. Ascent runs the extract afresh at every step, so it writes a file per step numbered by the Blueprint's `state/cycle`: the path is `$MIO_CATALYST_OUTPUT` when it holds `{step}`, else `ascent_{step}.vtu`. A multi-domain node writes its first domain, as `from_conduit` reads it.

## Reading it back

```python
for time, mesh in meshioplusplus.read_sequence("run.vtkhdf"):
    ...
```

## Checked against

The script is tested against Conduit nodes stood in by nested dicts in `tests/python/test_contrib_routes.py` (a VTKHDF series of three steps read back by `read_sequence`, the refused grid change, a channel named and fanned out, an Ascent cycle). The run inside ParaView Catalyst and inside Ascent is listed in the [roadmap](../roadmap.md#awaiting-a-licensed-run). The references are the [Catalyst 2 protocol](https://docs.paraview.org/en/latest/Catalyst/blueprints.html) and [Ascent's Python extract](https://ascent.readthedocs.io/en/latest/Actions/Extracts.html).
