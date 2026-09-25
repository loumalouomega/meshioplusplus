# MSC Marc `.t16`

Marc's binary post file is read through PyPost (`py_post`), the Python module Marc/Mentat ships. The route is `contrib/marc_t16/t16_to_vtu.py`. Where the formatted twin is available (`POST` with a `.t19`), meshio++ reads that natively ([Marc](../formats/marc.md)).

## Running it

Copy `t16_to_vtu.py` and `contrib/common/mio_vtu_series.py` into one directory, then, with the Python Mentat uses (its `bin` directory holds `py_post`):

```bash
python t16_to_vtu.py job.t16                     # job.pvd + job/
python t16_to_vtu.py job.t16 out.pvd --increments 1,10,20
```

## What it writes

- **Increments → steps.** Post file position 0 holds the model; every later position is a step at its time (`p.time`), made unique and ascending where two increments share one.
- **Node data.** The undeformed coordinates; every node vector (`Displacement`, `Reaction Force`…) and node scalar as point data.
- **Element data.** Every element scalar and tensor, averaged over the element's nodes, as cell data; a tensor has six components in the order 11, 22, 33, 12, 23, 13.
- **Ids and sets.** `marc:node_id` and `marc:element_id`, and every set as a `set:<name>` 0/1 mask.
- **Elements.** Types map as the native `.t19` reader maps them (Volume B; the extra nodes of Herrmann, bubble and generalized-plane-strain types dropped); others are skipped with a count. A remeshed increment carries its own mesh.

## Checked against

The script's logic is tested against a stand-in of PyPost (`tests/python/test_contrib_routes.py`); the run against Marc itself is listed in the [roadmap](../roadmap.md#awaiting-a-licensed-run). The reference is the Marc Python Reference's PyPost chapter.
