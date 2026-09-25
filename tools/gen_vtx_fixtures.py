"""Generate the DOLFINx VTX (``.bp``) fixtures under ``tests/python/meshes/vtx``.

Run inside a DOLFINx environment (DOLFINx 0.11 with ADIOS2 from conda-forge was
used; the versions are recorded in the fixtures' README)::

    python tools/gen_vtx_fixtures.py serial <outdir>
    mpirun -n 2 python tools/gen_vtx_fixtures.py parallel <outdir>
    python tools/gen_vtx_fixtures.py complex <outdir>   # plain adios2, no DOLFINx

``serial`` writes

- ``heat.bp`` -- the heat equation (P1 on triangles) stepped with backward Euler,
  three written steps, ``VTXMeshPolicy.reuse`` (the mesh is in step 0 only), BP5;
- ``elastic.bp`` -- a P2 vector field and a DG0 cell field on tetrahedra, two
  steps, ``VTXMeshPolicy.update`` (the mesh is rewritten every step), BP4;
- ``quads.bp`` -- a mesh-only file (no function) of a quadrilateral mesh, BP5.

``parallel`` writes ``heat_np2.bp``: the heat run on two ranks, so every step
holds two blocks with their own local point numbering and ghost points.

``complex`` writes ``complex.bp`` with the ``adios2`` package directly, in the
layout DOLFINx uses for a complex-valued function (``<name>_real`` and
``<name>_imag`` arrays); a complex DOLFINx build is a separate PETSc flavour.

``plain`` writes ``plain.bp``, an ADIOS2 file with no ``vtk.xml`` schema (the
reader must refuse it by name).

``raw`` records, for every ``<name>.bp``, the function arrays of every step as
the ``adios2`` package reads them (the rank blocks concatenated) in
``<name>.raw.npz`` (``k<k>/<name>``): the oracle for the steps ParaView cannot
read (see ``tools/gen_vtx_reference.py``).
"""

import sys
from pathlib import Path

import numpy as np


def _heat(comm, out, name, engine):
    import ufl
    from dolfinx import fem, mesh
    from dolfinx.fem.petsc import LinearProblem
    from dolfinx.io import VTXMeshPolicy, VTXWriter

    domain = mesh.create_unit_square(comm, 4, 4, mesh.CellType.triangle)
    V = fem.functionspace(domain, ("Lagrange", 1))
    u_n = fem.Function(V, name="u")
    u_n.interpolate(lambda x: np.exp(-20.0 * ((x[0] - 0.5) ** 2 + (x[1] - 0.5) ** 2)))
    dt = 0.01
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = u * v * ufl.dx + dt * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = u_n * v * ufl.dx
    problem = LinearProblem(
        a,
        L,
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        petsc_options_prefix="heat_",
    )
    with VTXWriter(
        comm, out / name, [u_n], engine=engine, mesh_policy=VTXMeshPolicy.reuse
    ) as vtx:
        t = 0.0
        vtx.write(t)
        for _ in range(2):
            t += dt
            uh = problem.solve()
            u_n.x.array[:] = uh.x.array
            vtx.write(t)


def serial(out):
    from dolfinx import fem, mesh
    from dolfinx.io import VTXMeshPolicy, VTXWriter
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    _heat(comm, out, "heat.bp", "BP5")

    domain = mesh.create_unit_cube(comm, 2, 2, 1, mesh.CellType.tetrahedron)
    V = fem.functionspace(domain, ("Lagrange", 2, (3,)))
    W = fem.functionspace(domain, ("DG", 0))
    disp = fem.Function(V, name="displacement")
    rho = fem.Function(W, name="density")
    with VTXWriter(
        comm,
        out / "elastic.bp",
        [disp, rho],
        engine="BP4",
        mesh_policy=VTXMeshPolicy.update,
    ) as vtx:
        for k, t in enumerate((0.0, 0.5)):
            disp.interpolate(
                lambda x, k=k: np.vstack((k * x[0] * x[1], -x[2], x[0] + k))
            )
            rho.interpolate(lambda x, k=k: 1.0 + k + x[0])
            vtx.write(t)

    quads = mesh.create_unit_square(comm, 2, 3, mesh.CellType.quadrilateral)
    with VTXWriter(comm, out / "quads.bp", quads, engine="BP5") as vtx:
        vtx.write(0.0)


def parallel(out):
    from mpi4py import MPI

    _heat(MPI.COMM_WORLD, out, "heat_np2.bp", "BP5")


def complex_(out):
    """A two-triangle mesh with one complex P1 field, written like DOLFINx does."""
    import adios2

    schema = (
        '<VTKFile type="UnstructuredGrid" version="0.1">\n'
        "  <UnstructuredGrid>\n"
        '    <Piece NumberOfPoints="NumberOfNodes" NumberOfCells="NumberOfEntities">\n'
        "      <Points>\n"
        '        <DataArray Name="geometry" />\n'
        "      </Points>\n"
        "      <Cells>\n"
        '        <DataArray Name="connectivity" />\n'
        '        <DataArray Name="types" />\n'
        "      </Cells>\n"
        "      <PointData>\n"
        '        <DataArray Name="vtkOriginalPointIds" />\n'
        '        <DataArray Name="vtkGhostType" />\n'
        '        <DataArray Name="p_real" />\n'
        '        <DataArray Name="p_imag" />\n'
        '        <DataArray Name="TIME">\n'
        "          step\n"
        "        </DataArray>\n"
        "      </PointData>\n"
        "    </Piece>\n"
        "  </UnstructuredGrid>\n"
        "</VTKFile>\n"
    )
    x = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=np.float64)
    conn = np.array([[3, 0, 1, 2], [3, 0, 2, 3]], dtype=np.int64)
    with adios2.Stream(str(out / "complex.bp"), "w") as fw:
        fw.write_attribute("vtk.xml", schema)
        for k, t in enumerate((0.0, 1.0)):
            fw.begin_step()
            fw.write("step", float(t))
            fw.write("geometry", x, [], [], list(x.shape))
            fw.write("NumberOfNodes", np.array([4], dtype=np.uint32), [], [], [1])
            fw.write("NumberOfCells", np.array([2], dtype=np.uint32), [], [], [1])
            fw.write("types", np.array([5], dtype=np.uint32), [], [], [1])
            fw.write("connectivity", conn, [], [], list(conn.shape))
            fw.write("vtkOriginalPointIds", np.arange(4, dtype=np.int64), [], [], [4])
            fw.write("vtkGhostType", np.zeros(4, dtype=np.uint8), [], [], [4])
            p = (k + 1) * (x[:, 0] + 1j * x[:, 1])
            fw.write("p_real", np.ascontiguousarray(p.real), [], [], [4])
            fw.write("p_imag", np.ascontiguousarray(p.imag), [], [], [4])
            fw.end_step()


def plain(out):
    import adios2

    with adios2.Stream(str(out / "plain.bp"), "w") as fw:
        fw.begin_step()
        fw.write("x", np.arange(4, dtype=np.float64), [], [], [4])
        fw.end_step()


_MESH_VARIABLES = {
    "step",
    "geometry",
    "connectivity",
    "types",
    "NumberOfNodes",
    "NumberOfEntities",
    "NumberOfCells",
    "vtkOriginalPointIds",
    "vtkGhostType",
}


def raw(out):
    import adios2

    for bp in sorted(out.glob("*.bp")):
        if bp.name == "plain.bp":
            continue
        arrays = {}
        with adios2.FileReader(str(bp)) as reader:
            for name in reader.available_variables():
                if name in _MESH_VARIABLES:
                    continue
                infos = reader.all_blocks_info(name)
                for k, blocks in enumerate(infos):
                    var = reader.inquire_variable(name)
                    parts = []
                    for index, info in enumerate(blocks):
                        var.set_step_selection([k, 1])
                        var.set_block_selection(index)
                        count = [int(c) for c in info["Count"].split(",")]
                        parts.append(np.asarray(reader.read(var)).reshape(count))
                    data = np.concatenate(parts)
                    if data.ndim == 2 and data.shape[1] == 1:
                        data = data[:, 0]
                    arrays[f"k{k}/{name}"] = data
        np.savez_compressed(bp.with_suffix(".raw.npz"), **arrays)


if __name__ == "__main__":
    mode, out = sys.argv[1], Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    {
        "serial": serial,
        "parallel": parallel,
        "complex": complex_,
        "plain": plain,
        "raw": raw,
    }[mode](out)
