# Femap `.modfem`

A Femap model file is Femap's private database, readable only by Femap. The route is to have Femap write a neutral file, which meshio++ reads natively with its groups and output sets ([Femap](../formats/femap.md)): *File → Export → Femap Neutral* by hand, or `contrib/femap/export_neutral.py` to script it.

## Running it

On Windows, with Femap installed and `pywin32` in any Python:

```bash
python export_neutral.py model.modfem                   # model.neu
python export_neutral.py model.modfem out.neu --version 2401
```

It attaches to a running Femap (else starts one) through the COM API, opens the model with `feFileOpen` and writes it with `feFileWriteNeutral`: the whole model (group 0), nodes, elements, properties, materials, groups and output sets, no geometry, in the neutral file version asked for (the default is Femap 2401's, which meshio++ reads).

`feFileWriteNeutral`'s argument list is the one to check against your Femap version's API help first: it was taken from the API reference, not from a run, and the run is listed in the [roadmap](../roadmap.md#awaiting-a-licensed-run). The script's own logic is tested against a stand-in (`tests/python/test_contrib_routes.py`).
