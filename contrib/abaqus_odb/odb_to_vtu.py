# MIT License -- part of meshio++ (https://github.com/loumalouomega/meshioplusplus).
"""Export an Abaqus output database (``.odb``) to a ParaView series meshio++ reads.

Run with Abaqus's own Python (the ODB API needs a licensed install)::

    abaqus python odb_to_vtu.py job.odb [job.pvd] [--fields U,S] [--steps Step-1]

``mio_vtu_series.py`` (from ``contrib/common``) must sit next to this script.
The output is ``job.pvd`` plus a ``job/`` directory of ``.vtu`` files:

- every part instance is a ``.pvd`` part named after the instance, so
  ``meshioplusplus.read("job.pvd")`` merges them with one region per instance;
- every frame is a time step: a time-domain frame at the step's total time at
  its start plus the frame value (frame 0 of a later step, the previous step's
  last state, is skipped); a frequency or modal frame (whose frame value is a
  frequency or an eigenvalue, not a time) one after the previous step's time, so
  the times keep ascending;
- field outputs by position: ``NODAL`` -> point data, ``INTEGRATION_POINT`` /
  ``CENTROID`` / ``WHOLE_ELEMENT`` -> cell data (Abaqus's own centroid subset,
  averaged over section points), ``ELEMENT_NODAL`` only -> point data averaged
  over the elements sharing a node;
- node and element labels as ``abaqus:node_label`` / ``abaqus:element_label``,
  and every node and element set (instance and assembly level) as a
  ``set:<name>`` 0/1 mask.

Runs on Python 2.7 (Abaqus 2023 and before) and 3 (Abaqus 2024 and later).
See ``doc/routes/abaqus_odb.md``.
"""

from __future__ import division, print_function

import argparse
import os
import sys

import numpy as np

# mio_vtu_series.py: copied next to this script, or in the repository's
# contrib/common.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_HERE, os.path.join(os.path.dirname(_HERE), "common")]

from mio_vtu_series import PvdSeries  # noqa: E402

# VTK cell types by node count, per element family (Abaqus lists the nodes of
# every shape below in VTK's order, as meshio++'s .inp and .fil readers read).
_SOLID = {4: 10, 6: 13, 8: 12, 10: 24, 15: 26, 20: 25}
_PLANAR = {3: 5, 4: 9, 6: 22, 8: 23, 9: 28}
_LINE = {2: 3, 3: 21}
_FAMILIES = (
    (("C3D", "DC3D", "AC3D", "COH3D", "SC"), _SOLID),
    (
        (
            "CPE",
            "CPS",
            "CAX",
            "CGAX",
            "DC2D",
            "DCAX",
            "COH2D",
            "STRI",
            "M3D",
            "R3D",
            "SFM3D",
            "S",
        ),
        _PLANAR,
    ),
    (("T2D", "T3D", "B2", "B3", "PIPE", "R2D", "RB2D", "RB3D", "DC1D"), _LINE),
)


def vtk_type(abaqus_type, nodes):
    """The VTK cell type of an Abaqus element, or None when it has no cell."""
    name = abaqus_type.upper()
    for prefixes, by_count in _FAMILIES:
        for prefix in prefixes:
            if not name.startswith(prefix):
                continue
            # A shell is S followed by a digit (S4R, S8R5); SC8R is a solid.
            if prefix == "S" and not (len(name) > 1 and name[1].isdigit()):
                continue
            return by_count.get(nodes)
    return None


class InstanceMesh(object):
    """One part instance: points, cell blocks and the label maps."""

    def __init__(self, instance):
        self.name = instance.name
        labels = [n.label for n in instance.nodes]
        self.node_index = dict((label, i) for i, label in enumerate(labels))
        points = np.zeros((len(labels), 3))
        for i, node in enumerate(instance.nodes):
            xyz = tuple(node.coordinates)
            points[i, : len(xyz)] = xyz
        self.points = points
        self.node_labels = np.array(labels, dtype=np.int64)

        by_type = {}
        order = []
        skipped = {}
        for element in instance.elements:
            conn = element.connectivity
            vtk = vtk_type(element.type, len(conn))
            if vtk is None:
                skipped[element.type] = skipped.get(element.type, 0) + 1
                continue
            if vtk not in by_type:
                by_type[vtk] = []
                order.append(vtk)
            by_type[vtk].append(
                (element.label, [self.node_index[label] for label in conn])
            )
        for etype, count in sorted(skipped.items()):
            print(
                "odb_to_vtu: instance %s: %d %s elements have no VTK cell, skipped"
                % (self.name, count, etype)
            )
        self.blocks = []
        element_labels = []
        for vtk in order:
            rows = by_type[vtk]
            self.blocks.append((vtk, np.array([r[1] for r in rows], dtype=np.int64)))
            element_labels.extend(r[0] for r in rows)
        self.element_labels = np.array(element_labels, dtype=np.int64)
        self.element_index = dict((label, i) for i, label in enumerate(element_labels))

    def node_mask(self, labels):
        mask = np.zeros(len(self.points), dtype=np.uint8)
        for label in labels:
            if label in self.node_index:
                mask[self.node_index[label]] = 1
        return mask

    def element_mask(self, labels):
        mask = np.zeros(len(self.element_labels), dtype=np.uint8)
        for label in labels:
            if label in self.element_index:
                mask[self.element_index[label]] = 1
        return mask


def _set_members(odb_set, instance_name, attribute):
    """Labels of an instance's members of a node or element set. An assembly
    set holds one sequence per instance (``instanceNames`` gives the order);
    an instance set holds the members directly."""
    members = getattr(odb_set, attribute)
    names = getattr(odb_set, "instanceNames", None)
    if names and len(members) and not hasattr(members[0], "label"):
        out = []
        for name, group in zip(names, members):
            if name == instance_name:
                out.extend(m.label for m in group)
        return out
    return [
        m.label
        for m in members
        if getattr(m, "instanceName", None) in (None, "", instance_name)
    ]


def set_masks(odb, instance, mesh):
    point_sets, cell_sets = {}, {}
    sources = [instance]
    if hasattr(odb, "rootAssembly"):
        sources.append(odb.rootAssembly)
    for source in sources:
        for name in sorted(getattr(source, "nodeSets", {}).keys()):
            labels = _set_members(source.nodeSets[name], mesh.name, "nodes")
            if labels:
                point_sets["set:" + name] = mesh.node_mask(labels)
        for name in sorted(getattr(source, "elementSets", {}).keys()):
            labels = _set_members(source.elementSets[name], mesh.name, "elements")
            if labels:
                cell_sets["set:" + name] = mesh.element_mask(labels)
    return point_sets, cell_sets


def _entries(subset):
    """(labels, data rows, is_nodal) of a field subset, through bulkDataBlocks
    where the ODB API has it, else value by value."""
    blocks = getattr(subset, "bulkDataBlocks", None)
    if blocks:
        for block in blocks:
            data = np.asarray(block.data, dtype=np.float64)
            node_labels = getattr(block, "nodeLabels", None)
            element_labels = getattr(block, "elementLabels", None)
            if node_labels is not None and len(node_labels):
                yield np.asarray(node_labels), data, True
            else:
                yield np.asarray(element_labels), data, False
        return
    labels, rows, nodal = [], [], None
    for value in subset.values:
        node_label = getattr(value, "nodeLabel", None)
        if node_label is not None:
            labels.append(node_label)
            nodal = True
        else:
            labels.append(value.elementLabel)
            nodal = False
        rows.append(np.atleast_1d(np.asarray(value.data, dtype=np.float64)))
    if labels:
        yield np.asarray(labels), np.asarray(rows), nodal


def _accumulate(n, index, entries):
    total, count, ncomp = None, np.zeros(n), None
    for labels, data, _ in entries:
        if data.ndim == 1:
            data = data[:, None]
        if total is None:
            ncomp = data.shape[1]
            total = np.zeros((n, ncomp))
        for label, row in zip(labels, data):
            i = index.get(int(label))
            if i is not None:
                total[i] += row
                count[i] += 1
    if total is None:
        return None
    with np.errstate(invalid="ignore", divide="ignore"):
        out = total / count[:, None]
    out[count == 0] = np.nan
    return out[:, 0] if ncomp == 1 else out


def field_arrays(field, instance, mesh, constants):
    """(point data, cell data) of one field output on one instance."""
    positions = set(loc.position for loc in field.locations)
    point, cell = {}, {}
    if constants.NODAL in positions:
        sub = field.getSubset(region=instance, position=constants.NODAL)
        values = _accumulate(len(mesh.points), mesh.node_index, _entries(sub))
        if values is not None:
            point[field.name] = values
    elif positions & set(
        (constants.INTEGRATION_POINT, constants.CENTROID, constants.WHOLE_ELEMENT)
    ):
        position = (
            constants.WHOLE_ELEMENT
            if constants.WHOLE_ELEMENT in positions
            else constants.CENTROID
        )
        sub = field.getSubset(region=instance, position=position)
        values = _accumulate(
            len(mesh.element_labels), mesh.element_index, _entries(sub)
        )
        if values is not None:
            cell[field.name] = values
    elif constants.ELEMENT_NODAL in positions:
        sub = field.getSubset(region=instance, position=constants.ELEMENT_NODAL)
        values = _accumulate(len(mesh.points), mesh.node_index, _entries(sub))
        if values is not None:
            point[field.name] = values
    return point, cell


def frames(odb, step_names):
    """(time, step name, frame) for every exported frame, times ascending."""
    last = None
    for step_name in odb.steps.keys():
        if step_names and step_name not in step_names:
            continue
        step = odb.steps[step_name]
        for k, frame in enumerate(step.frames):
            domain = str(getattr(frame, "domain", "TIME"))
            if domain.endswith("TIME"):
                t = float(step.totalTime) + float(frame.frameValue)
                if k == 0 and last is not None and t <= last:
                    continue  # the previous step's last state
            else:
                t = (last if last is not None else -1.0) + 1.0
            if last is not None and t <= last:
                t = last + 1.0
            last = t
            yield t, step_name, frame


def export(odb, out, fields=None, steps=None, constants=None):
    if constants is None:
        import abaqusConstants as constants
    instances = odb.rootAssembly.instances
    meshes = []
    for name in instances.keys():
        mesh = InstanceMesh(instances[name])
        if len(mesh.points) and mesh.blocks:
            meshes.append(
                (instances[name], mesh, set_masks(odb, instances[name], mesh))
            )
    if not meshes:
        raise SystemExit("odb_to_vtu: no instance has a cell meshio++ can hold")
    # Every set on every instance (0 where it has no member), so the merged
    # parts carry a 0/1 mask rather than NaN.
    point_names = set().union(*(ps for _, _, (ps, _) in meshes))
    cell_names = set().union(*(cs for _, _, (_, cs) in meshes))
    for _, mesh, (point_sets, cell_sets) in meshes:
        for name in point_names - set(point_sets):
            point_sets[name] = np.zeros(len(mesh.points), dtype=np.uint8)
        for name in cell_names - set(cell_sets):
            cell_sets[name] = np.zeros(len(mesh.element_labels), dtype=np.uint8)
    written = 0
    with PvdSeries(out) as series:
        for t, _, frame in frames(odb, steps):
            for part, (instance, mesh, (point_sets, cell_sets)) in enumerate(meshes):
                point = {"abaqus:node_label": mesh.node_labels}
                cell = {"abaqus:element_label": mesh.element_labels}
                point.update(point_sets)
                cell.update(cell_sets)
                for name in frame.fieldOutputs.keys():
                    if fields and name not in fields:
                        continue
                    p, c = field_arrays(
                        frame.fieldOutputs[name], instance, mesh, constants
                    )
                    point.update(p)
                    cell.update(c)
                series.add(
                    t,
                    mesh.points,
                    mesh.blocks,
                    point_data=point,
                    cell_data=cell,
                    part=part,
                    name=mesh.name,
                )
            written += 1
    return written


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("odb")
    parser.add_argument("out", nargs="?", help="the .pvd to write (default: <odb>.pvd)")
    parser.add_argument("--fields", help="comma-separated field outputs (default: all)")
    parser.add_argument("--steps", help="comma-separated step names (default: all)")
    args = parser.parse_args(argv)
    out = args.out or os.path.splitext(args.odb)[0] + ".pvd"
    fields = set(args.fields.split(",")) if args.fields else None
    steps = set(args.steps.split(",")) if args.steps else None

    from odbAccess import openOdb

    odb = openOdb(args.odb, readOnly=True)
    try:
        count = export(odb, out, fields, steps)
    finally:
        odb.close()
    print("odb_to_vtu: wrote %d frames to %s" % (count, out))


if __name__ == "__main__":
    main()
