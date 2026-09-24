// The parallel MFEM fixtures of tests/python/meshes/mfem/parallel, built and run
// by tools/gen_mfem_parallel_fixtures.py (not part of meshio++; needs MFEM built
// with MPI and hypre, both BSD/MIT-style licensed).
//
//   mfem_parallel_driver <serial mesh> <out prefix> <curve order (0 = none)> <u order>
//
// Refines the mesh once, partitions it by element centroid x over the MPI ranks
// and writes, per rank: <out>.mesh.%06d (ParMesh::Save, the files GLVis reads),
// <out>.pmesh.%06d (ParMesh::ParPrint, with the communication groups),
// <out>.u.%06d (an H1 field) and MFEM's own high-order VTU (ParaViewDataCollection)
// under <out>_ref/.
#include "mfem.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace mfem;

static void warp(const Vector& x, Vector& y) {
    y = x;
    const int d = x.Size();
    y(0) = x(0) + 0.1 * std::sin(1.3 * (d > 1 ? x(1) : x(0)));
    if (d > 1)
        y(1) = x(1) + 0.07 * std::cos(0.9 * x(0));
    if (d > 2)
        y(2) = x(2) + 0.05 * std::sin(x(0) + x(1));
}

static double field(const Vector& x) {
    double s = 0;
    for (int i = 0; i < x.Size(); ++i)
        s += x(i);
    return std::sin(x(0)) + std::pow(x(x.Size() - 1), 3) + 0.5 * std::cos(s);
}

int main(int argc, char* argv[]) {
    Mpi::Init(argc, argv);
    const int np = Mpi::WorldSize(), rank = Mpi::WorldRank();
    const char* src = argv[1];
    const std::string out = argv[2];
    const int curve = std::atoi(argv[3]);
    const int u_order = std::atoi(argv[4]);
    Mesh serial(src, 1, 1);
    serial.UniformRefinement();
    // Partition by the element centroid's x.
    Array<int> part(serial.GetNE());
    double lo = 1e300, hi = -1e300;
    std::vector<double> cx(serial.GetNE());
    for (int e = 0; e < serial.GetNE(); ++e) {
        Vector c;
        serial.GetElementCenter(e, c);
        cx[e] = c(0);
        lo = std::min(lo, c(0));
        hi = std::max(hi, c(0));
    }
    for (int e = 0; e < serial.GetNE(); ++e)
        part[e] = std::min(np - 1, int((cx[e] - lo) / (hi - lo + 1e-12) * np));
    ParMesh pmesh(MPI_COMM_WORLD, serial, part.GetData());
    if (curve > 0) {
        pmesh.SetCurvature(curve, false, -1, Ordering::byVDIM);
        VectorFunctionCoefficient w(pmesh.SpaceDimension(), warp);
        // Project into a copy: projecting onto the nodes in place would read
        // nodes already moved, in each rank's own element order.
        ParGridFunction* nodes = dynamic_cast<ParGridFunction*>(pmesh.GetNodes());
        ParGridFunction warped(nodes->ParFESpace());
        warped.ProjectCoefficient(w);
        *nodes = warped;
    }
    H1_FECollection fec(u_order, pmesh.Dimension());
    ParFiniteElementSpace fes(&pmesh, &fec);
    ParGridFunction u(&fes);
    FunctionCoefficient f(field);
    u.ProjectCoefficient(f);
    pmesh.Save(out + ".mesh", 17);
    {
        std::ostringstream name;
        name << out << ".pmesh." << std::setfill('0') << std::setw(6) << rank;
        std::ofstream ofs(name.str());
        ofs.precision(17);
        pmesh.ParPrint(ofs);
    }
    u.Save((out + ".u").c_str(), 17);
    const int order = std::max(u_order, curve > 0 ? curve : 1);
    ParaViewDataCollection dc("ref", &pmesh);
    dc.SetPrefixPath(out + "_ref");
    dc.SetLevelsOfDetail(order);
    dc.SetHighOrderOutput(true);
    dc.SetDataFormat(VTKFormat::ASCII);
    dc.SetPrecision(17);
    dc.RegisterField("u", &u);
    dc.Save();
    if (rank == 0)
        std::cout << out << ": " << np << " ranks, " << serial.GetNE() << " elements, "
                  << serial.GetNV() << " vertices, order " << order << "\n";
    return 0;
}
