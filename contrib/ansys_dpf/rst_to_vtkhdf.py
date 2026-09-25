# MIT License -- part of meshio++ (https://github.com/loumalouomega/meshioplusplus).
"""Export Ansys results through PyDPF to a VTKHDF time series meshio++ reads.

For the result files and quantities meshio++'s native reader does not cover
(``doc/formats/ansys_rst.md``: ``.rth``/``.rmg``, compressed ``/FCOMP``
records, results DPF derives), let Ansys's own Data Processing Framework read
them. Runs in an ordinary Python with ``ansys-dpf-core``, ``pyvista`` and
meshio++ installed, next to a DPF Server (shipped with Ansys 2021 R1 and
later, or the standalone DPF Server)::

    python rst_to_vtkhdf.py file.rst [file.vtkhdf] [--results displacement,stress]

Every result set is a time step at its time or frequency. Nodal results are
point data; elemental ones cell data; elemental-nodal ones (stress, strain) are
averaged to the nodes by DPF, as the Mechanical and MAPDL plots show them.
Node and element numbers ride along as ``dpf:node_id`` / ``dpf:element_id``.
See ``doc/routes/ansys_dpf.md``.
"""

import argparse
import os

import numpy as np

import meshioplusplus

DEFAULT_RESULTS = ("displacement", "stress", "elastic_strain", "temperature")


def _available(model):
    return {r.name: r for r in model.metadata.result_info.available_results}


def result_fields(model, dpf, name, set_ids):
    """``{set id: (location, ids, data)}`` for one result, over all sets."""
    info = _available(model)[name]
    op = model.operator(info.operator_name)
    op.connect(0, dpf.Scoping(ids=list(set_ids), location=dpf.locations.time_freq))
    if info.native_location == dpf.locations.elemental_nodal:
        op.connect(9, dpf.locations.nodal)
    fc = op.get_output(0, dpf.types.fields_container)
    out = {}
    for set_id in set_ids:
        fields = fc.get_fields({"time": set_id})
        if not fields:
            continue
        location = fields[0].location
        ids = np.concatenate([np.asarray(f.scoping.ids) for f in fields])
        data = np.concatenate([np.asarray(f.data, dtype=np.float64) for f in fields])
        out[set_id] = (location, ids, data)
    return out


def _scatter(n, id_to_index, ids, data):
    shape = (n,) + tuple(data.shape[1:])
    out = np.full(shape, np.nan)
    for i, entity in enumerate(ids):
        k = id_to_index.get(int(entity))
        if k is not None:
            out[k] = data[i]
    return out


def steps(model, dpf, results):
    mesh = model.metadata.meshed_region
    grid = mesh.grid
    node_ids = np.asarray(mesh.nodes.scoping.ids)
    element_ids = np.asarray(mesh.elements.scoping.ids)
    node_index = {int(i): k for k, i in enumerate(node_ids)}
    element_index = {int(i): k for k, i in enumerate(element_ids)}
    support = model.metadata.time_freq_support
    times = np.asarray(support.time_frequencies.data, dtype=np.float64)
    set_ids = list(range(1, len(times) + 1))
    available = _available(model)
    fields = {
        name: result_fields(model, dpf, name, set_ids)
        for name in results
        if name in available
    }
    for k, set_id in enumerate(set_ids):
        g = grid.copy()
        g.point_data["dpf:node_id"] = node_ids
        g.cell_data["dpf:element_id"] = element_ids
        for name, by_set in fields.items():
            if set_id not in by_set:
                continue
            location, ids, data = by_set[set_id]
            if location == dpf.locations.nodal:
                g.point_data[name] = _scatter(len(node_ids), node_index, ids, data)
            elif location == dpf.locations.elemental:
                g.cell_data[name] = _scatter(len(element_ids), element_index, ids, data)
        yield float(times[k]), meshioplusplus.from_pyvista(g)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("rst")
    parser.add_argument(
        "out", nargs="?", help="the .vtkhdf to write (default: <rst>.vtkhdf)"
    )
    parser.add_argument(
        "--results",
        default=",".join(DEFAULT_RESULTS),
        help="comma-separated DPF result names (default: %(default)s)",
    )
    args = parser.parse_args(argv)
    out = args.out or os.path.splitext(args.rst)[0] + ".vtkhdf"

    from ansys.dpf import core as dpf

    model = dpf.Model(args.rst)
    results = [r for r in args.results.split(",") if r]
    written = meshioplusplus.write_sequence(out, steps(model, dpf, results))
    print("rst_to_vtkhdf: wrote %s" % written[0])


if __name__ == "__main__":
    main()
