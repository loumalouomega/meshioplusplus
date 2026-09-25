# Ansys results through PyDPF

meshio++ reads MAPDL `.rst`/`.rth` natively ([Ansys results](../formats/ansys_rst.md)). For what that reader does not cover — compressed `/FCOMP` records, the thermal-electromagnetic `.rmg`, the results Mechanical derives, Ansys formats with no MAPDL twin — Ansys's own Data Processing Framework reads the file, and `contrib/ansys_dpf/rst_to_vtkhdf.py` hands the result to meshio++.

## Running it

In an ordinary Python with `ansys-dpf-core`, `pyvista` and `meshioplusplus[h5py]` installed, and a DPF Server (part of Ansys 2021 R1 and later, or the standalone DPF Server):

```bash
python rst_to_vtkhdf.py file.rst                            # file.vtkhdf
python rst_to_vtkhdf.py file.rst out.vtkhdf --results displacement,stress,temperature
```

## What it writes

- **Result sets → steps**, at each set's time or frequency, as one [VTKHDF](../formats/vtkhdf.md) time series.
- **Results**: nodal results are point data and elemental ones cell data; elemental-nodal results (stress, strain) are averaged to the nodes by DPF itself, as MAPDL and Mechanical plot them. The defaults are `displacement`, `stress`, `elastic_strain` and `temperature`, each where the file has it; `--results` takes any DPF result name.
- **Ids**: `dpf:node_id` and `dpf:element_id`.
- The mesh is DPF's (`meshed_region.grid`), converted with `meshioplusplus.from_pyvista`.

## Checked against

The script's logic is tested against a stand-in of the DPF API (`tests/python/test_contrib_routes.py`); the run against a DPF Server is listed in the [roadmap](../roadmap.md#awaiting-a-licensed-run). The reference is the [PyDPF-Core documentation](https://dpf.docs.pyansys.com/).
