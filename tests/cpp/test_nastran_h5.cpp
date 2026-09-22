//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

/**
 * @file test_nastran_h5.cpp
 * @brief MSC Nastran HDF5 reader: the compound-table helpers, the element
 *        permutations, steps, and the error paths.
 *
 * The files are built here with the HDF5 C API in the MSC schema rather than read
 * from tests/python/meshes/nastran_h5/: this suite has no test-data path. The
 * Python suite runs the same reader over real MSC Nastran output.
 */

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <hdf5.h>

// Project includes
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/nastran_h5.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

namespace h5 = meshioplusplus::h5;
using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;

/** One member of a compound table: int64 or float64, `mWidth` wide (an ARRAY if > 1). */
struct Field {
    std::string mName;
    bool mReal;
    std::size_t mWidth;
    std::vector<double> mValues;  // rows * width, row-major
};

/** Writes a rank-1 compound dataset at `rPath`, creating the groups on the way. */
void write_table(hid_t file, const std::string& rPath, std::size_t Rows,
                 const std::vector<Field>& rFields) {
    std::vector<std::size_t> offsets;
    std::size_t size = 0;
    for (const Field& f : rFields) {
        offsets.push_back(size);
        size += 8 * f.mWidth;
    }
    h5::Hid type(H5Tcreate(H5T_COMPOUND, size), H5Tclose);
    for (std::size_t k = 0; k < rFields.size(); ++k) {
        const Field& f = rFields[k];
        const hid_t base = f.mReal ? H5T_NATIVE_DOUBLE : H5T_NATIVE_INT64;
        if (f.mWidth == 1) {
            H5Tinsert(type, f.mName.c_str(), offsets[k], base);
        } else {
            const hsize_t dim = f.mWidth;
            h5::Hid arr(H5Tarray_create2(base, 1, &dim), H5Tclose);
            H5Tinsert(type, f.mName.c_str(), offsets[k], arr);
        }
    }
    std::vector<unsigned char> buf(size * Rows);
    for (std::size_t r = 0; r < Rows; ++r)
        for (std::size_t k = 0; k < rFields.size(); ++k)
            for (std::size_t w = 0; w < rFields[k].mWidth; ++w) {
                unsigned char* dst = buf.data() + r * size + offsets[k] + 8 * w;
                const double v = rFields[k].mValues[r * rFields[k].mWidth + w];
                if (rFields[k].mReal) {
                    std::memcpy(dst, &v, 8);
                } else {
                    const std::int64_t i = static_cast<std::int64_t>(v);
                    std::memcpy(dst, &i, 8);
                }
            }
    const hsize_t n = Rows;
    h5::Hid space(H5Screate_simple(1, &n, nullptr), H5Sclose);
    h5::Hid lcpl(H5Pcreate(H5P_LINK_CREATE), H5Pclose);
    H5Pset_create_intermediate_group(lcpl, 1);
    h5::Hid d(H5Dcreate2(file, rPath.c_str(), type, space, lcpl, H5P_DEFAULT, H5P_DEFAULT),
              H5Dclose);
    ASSERT_TRUE(d.Valid()) << rPath;
    H5Dwrite(d, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
}

void write_version(hid_t file, const std::string& rVersion) {
    h5::Hid g(H5Gopen2(file, "/NASTRAN", H5P_DEFAULT), H5Gclose);
    h5::Hid t(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(t, rVersion.size());
    h5::Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    h5::Hid a(H5Acreate2(g, "VERSION", t, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    H5Awrite(a, t, rVersion.data());
}

// The unit cube's corners, GRID 1-8, then the twelve edge midpoints, GRID 101-112,
// in Nastran's CHEXA20 order: bottom ring, vertical edges, top ring.
const double kCorner[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                              {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
const int kNastranEdges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5},
                                  {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}};

void write_grids(hid_t file, bool WithMids) {
    std::vector<double> id;
    std::vector<double> x;
    for (int i = 0; i < 8; ++i) {
        id.push_back(i + 1);
        x.insert(x.end(), kCorner[i], kCorner[i] + 3);
    }
    if (WithMids)
        for (int e = 0; e < 12; ++e) {
            id.push_back(101 + e);
            for (int c = 0; c < 3; ++c)
                x.push_back(0.5 *
                            (kCorner[kNastranEdges[e][0]][c] + kCorner[kNastranEdges[e][1]][c]));
        }
    const std::size_t n = id.size();
    const std::vector<double> zero(n, 0.0);
    write_table(file, "/NASTRAN/INPUT/NODE/GRID", n,
                {{"ID", false, 1, id},
                 {"CP", false, 1, zero},
                 {"X", true, 3, x},
                 {"CD", false, 1, zero},
                 {"DOMAIN_ID", false, 1, std::vector<double>(n, 1.0)}});
}

/** A linear CHEXA (EID 7, PID 3) with DISPLACEMENT in two domains; GRID 3 absent from the second.
 */
std::string hex_with_results(const std::string& rVersion = "msc20200") {
    const std::string path = mt::temp_path(".h5");
    h5::Hid file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
    write_grids(file, false);
    write_version(file, rVersion);
    std::vector<double> g{1, 2, 3, 4, 5, 6, 7, 8};
    g.resize(20, 0.0);
    write_table(file, "/NASTRAN/INPUT/ELEMENT/CHEXA", 1,
                {{"EID", false, 1, {7}},
                 {"PID", false, 1, {3}},
                 {"G", false, 20, g},
                 {"DOMAIN_ID", false, 1, {1}}});
    write_table(file, "/NASTRAN/INPUT/PROPERTY/PSOLID", 1,
                {{"PID", false, 1, {3}}, {"MID", false, 1, {1}}, {"DOMAIN_ID", false, 1, {1}}});
    write_table(file, "/NASTRAN/RESULT/DOMAINS", 2,
                {{"ID", false, 1, {1, 2}},
                 {"SUBCASE", false, 1, {1, 1}},
                 {"ANALYSIS", false, 1, {6, 6}},
                 {"TIME_FREQ_EIGR", true, 1, {0.5, 1.5}},
                 {"MODE", false, 1, {0, 0}}});
    std::vector<double> id;
    std::vector<double> x;
    std::vector<double> rz;
    std::vector<double> dom;
    for (int d = 1; d <= 2; ++d)
        for (int n = 1; n <= 8; ++n) {
            if (d == 2 && n == 3)
                continue;
            id.push_back(n);
            x.push_back(d == 1 ? n : -n);
            rz.push_back(d - 1);
            dom.push_back(d);
        }
    const std::size_t rows = id.size();
    const std::vector<double> zero(rows, 0.0);
    write_table(file, "/NASTRAN/RESULT/NODAL/DISPLACEMENT", rows,
                {{"ID", false, 1, id},
                 {"X", true, 1, x},
                 {"Y", true, 1, zero},
                 {"Z", true, 1, zero},
                 {"RX", true, 1, zero},
                 {"RY", true, 1, zero},
                 {"RZ", true, 1, rz},
                 {"DOMAIN_ID", false, 1, dom}});
    write_table(file, "/INDEX/NASTRAN/RESULT/NODAL/DISPLACEMENT", 2,
                {{"DOMAIN_ID", false, 1, {1, 2}},
                 {"POSITION", false, 1, {0, 8}},
                 {"LENGTH", false, 1, {8, 7}}});
    return path;
}

}  // namespace

TEST(NastranH5, CompoundHelpersReadOneMemberOverARowRange) {
    const std::string path = hex_with_results();
    h5::Hid file = h5::open_file_read(path);
    const auto members = h5::compound_members(file, "/NASTRAN/INPUT/NODE/GRID");
    ASSERT_EQ(members.size(), 5u);
    EXPECT_EQ(members[2].mName, "X");
    EXPECT_TRUE(members[2].mNumeric);
    EXPECT_EQ(members[2].mDtype, DType::Float64);
    EXPECT_EQ(members[2].mDims, std::vector<std::size_t>{3});
    EXPECT_TRUE(members[0].mDims.empty());

    const auto x = h5::read_compound_member(file, "/NASTRAN/INPUT/NODE/GRID", "X", 1, 2);
    ASSERT_EQ(x.Shape(), (std::vector<std::size_t>{2, 3}));
    EXPECT_DOUBLE_EQ(x.As<double>()[0], 1.0);  // GRID 2 = (1, 0, 0)
    EXPECT_DOUBLE_EQ(x.As<double>()[4], 1.0);  // GRID 3 = (1, 1, 0)

    const DType as = DType::Float64;
    const auto id = h5::read_compound_member(file, "/NASTRAN/INPUT/NODE/GRID", "ID", 0, 8, &as);
    EXPECT_EQ(id.Dtype(), DType::Float64);
    EXPECT_DOUBLE_EQ(id.As<double>()[7], 8.0);

    EXPECT_THROW(h5::read_compound_member(file, "/NASTRAN/INPUT/NODE/GRID", "NOPE", 0, 1),
                 ReadError);
    EXPECT_THROW(h5::read_compound_member(file, "/NASTRAN/INPUT/NODE/GRID", "X", 7, 2), ReadError);
    std::filesystem::remove(path);
}

TEST(NastranH5, Hexa20MidNodesSitOnTheirEdges) {
    const std::string path = mt::temp_path(".h5");
    {
        h5::Hid file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
        write_grids(file, true);
        std::vector<double> g{1, 2, 3, 4, 5, 6, 7, 8};
        for (int e = 0; e < 12; ++e)
            g.push_back(101 + e);
        write_table(file, "/NASTRAN/INPUT/ELEMENT/CHEXA", 1,
                    {{"EID", false, 1, {1}}, {"PID", false, 1, {1}}, {"G", false, 20, g}});
    }
    const Mesh mesh = meshioplusplus::read_nastran_h5(path);
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "hexahedron20");
    // meshio++ (VTK) order: bottom ring, top ring, vertical edges.
    const int vtk_edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                  {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    const auto* conn = mesh.Cells(0).Conn().As<std::int64_t>();
    const double* p = mesh.Points().As<double>();
    for (int e = 0; e < 12; ++e)
        for (int c = 0; c < 3; ++c)
            EXPECT_DOUBLE_EQ(p[conn[8 + e] * 3 + c], 0.5 * (p[conn[vtk_edges[e][0]] * 3 + c] +
                                                            p[conn[vtk_edges[e][1]] * 3 + c]))
                << "edge " << e;
    std::filesystem::remove(path);
}

TEST(NastranH5, EachDomainIsAStepWithNaNWhereATableHasNoRow) {
    const std::string path = hex_with_results();
    const Mesh first = meshioplusplus::read_nastran_h5(path);
    EXPECT_EQ(first.Cells(0).Type(), "hexahedron");
    EXPECT_DOUBLE_EQ(first.FieldData("meshio:time").As<double>()[0], 0.5);
    EXPECT_EQ(first.CellData("nastran:eid", 0).As<std::int64_t>()[0], 7);
    ASSERT_EQ(first.NumRegions(), 1u);
    EXPECT_EQ(first.Region(0).mName, "PSOLID_3");
    EXPECT_EQ(first.Region(0).mTag, 3);

    ReadOptions opts;
    opts.mTimeStep = -1;
    const Mesh second = meshioplusplus::read_nastran_h5(path, opts);
    EXPECT_EQ(second.FieldData("nastran:domain").As<std::int64_t>()[0], 2);
    const double* d = second.PointData("DISPLACEMENT").As<double>();
    EXPECT_TRUE(std::isnan(d[2 * 3]));  // GRID 3 has no row in domain 2
    EXPECT_DOUBLE_EQ(d[3 * 3], -4.0);
    EXPECT_DOUBLE_EQ(second.PointData("DISPLACEMENT_ROT").As<double>()[2], 1.0);

    opts.mTimeStep = 2;
    EXPECT_THROW(meshioplusplus::read_nastran_h5(path, opts), ReadError);

    const auto meta = meshioplusplus::read_nastran_h5_metadata(path);
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{0.5, 1.5}));
    EXPECT_EQ(meshioplusplus::sequence_num_steps(path, ""), 2u);
    std::filesystem::remove(path);
}

TEST(NastranH5, ArraysNarrowTheResultTables) {
    const std::string path = hex_with_results();
    ReadOptions opts;
    opts.mDataArrays = std::vector<std::string>{"DISPLACEMENT_ROT"};
    const Mesh mesh = meshioplusplus::read_nastran_h5(path, opts);
    EXPECT_EQ(mesh.NumPointData(), 1u);
    EXPECT_TRUE(mesh.HasPointData("DISPLACEMENT_ROT"));
    std::filesystem::remove(path);
}

TEST(NastranH5, OtherVendorsAndOtherFilesAreRefused) {
    const std::string nx = hex_with_results("NX Nastran 2019");
    EXPECT_THROW(meshioplusplus::read_nastran_h5(nx), ReadError);
    std::filesystem::remove(nx);

    const std::string other = mt::temp_path(".h5");
    {
        h5::Hid file(H5Fcreate(other.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
    }
    EXPECT_THROW(meshioplusplus::read_nastran_h5(other), ReadError);
    std::filesystem::remove(other);

    const std::string text = mt::temp_path(".h5");
    std::ofstream(text) << "hello\n";
    EXPECT_THROW(meshioplusplus::read_nastran_h5(text), ReadError);
    std::filesystem::remove(text);
}

TEST(NastranH5, AnUndefinedGridIsAnError) {
    const std::string path = mt::temp_path(".h5");
    {
        h5::Hid file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
        write_grids(file, false);
        write_table(file, "/NASTRAN/INPUT/ELEMENT/CROD", 1,
                    {{"EID", false, 1, {1}}, {"PID", false, 1, {1}}, {"G", false, 2, {1, 99}}});
    }
    EXPECT_THROW(meshioplusplus::read_nastran_h5(path), ReadError);
    std::filesystem::remove(path);
}

TEST(NastranH5, RegistryReachesItByExtension) {
    EXPECT_EQ(meshioplusplus::resolve_format("job.h5", ""), "nastran_h5");
    EXPECT_EQ(meshioplusplus::resolve_format("job.post.h5", ""), "gid");
    const std::string path = hex_with_results();
    const Mesh mesh = meshioplusplus::registry_read(path, meshioplusplus::resolve_format(path, ""),
                                                    ReadOptions{});
    EXPECT_TRUE(mesh.HasPointData("DISPLACEMENT"));
    std::filesystem::remove(path);
}

#endif  // MESHIOPLUSPLUS_HAS_HDF5
