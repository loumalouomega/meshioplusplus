// Generate the Tecplot SZL fixtures under tests/python/meshes/szplt.
//
// Every data set is described once and written twice: as .szplt through
// TecIO's writer (TecIO's new API writes only SZL), and as an ASCII Tecplot
// .dat twin written here, which meshio++ reads natively. The tests check that
// the TecIO-backed .szplt reader and the native ASCII reader agree.
//
// Build against a TecIO (the source Tecplot distributes, or a Tecplot 360
// install) and run in the output directory:
//
//   g++ -std=c++17 -I<tecio>/include tools/gen_szplt_fixtures.cpp
//       <tecio>/lib/libtecio.a -lpthread -o gen_szplt_fixtures
//   (cd tests/python/meshes/szplt && <path>/gen_szplt_fixtures)

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TECIO.h"

namespace {

void check(int Status, const char* pWhat) {
    if (Status != 0) {
        std::fprintf(stderr, "TecIO failed in %s\n", pWhat);
        std::exit(1);
    }
}

// TecIO's zone types, the ASCII ZONETYPE names.
const char* const kZoneTypes[] = {"ORDERED",         "FELINESEG",     "FETRIANGLE",
                                  "FEQUADRILATERAL", "FETETRAHEDRON", "FEBRICK"};

struct Zone {
    std::string mTitle;
    int32_t mType = 0;               // 0 ordered, else FE
    int64_t mI = 1, mJ = 1, mK = 1;  // ordered: I, J, K; FE: nodes, cells
    std::vector<int32_t> mLocation;  // per variable: 1 nodal, 0 cell-centred
    std::vector<int32_t> mShare;     // per variable: 1-based zone, or 0
    std::vector<int32_t> mPassive;   // per variable
    int32_t mConnShare = 0;
    bool mUnsteady = false;
    double mTime = 0.0;
    int32_t mStrand = 0;
    std::vector<std::vector<double>> mValues;  // per variable; empty when not owned
    std::vector<int32_t> mNodes;               // 1-based node map, when owned
};

struct DataSet {
    std::vector<std::string> mVariables;
    int32_t mVarType = 2;  // 1 float, 2 double
    std::vector<Zone> mZones;
};

void write_szplt(const std::string& rName, const DataSet& rD) {
    std::string vars;
    for (std::size_t v = 0; v < rD.mVariables.size(); ++v)
        vars += (v ? "," : "") + rD.mVariables[v];
    void* h = nullptr;
    check(tecFileWriterOpen(rName.c_str(), "meshio++ fixture", vars.c_str(), 1, 0, rD.mVarType,
                            nullptr, &h),
          "tecFileWriterOpen");
    const std::size_t nv = rD.mVariables.size();
    for (const Zone& rZ : rD.mZones) {
        std::vector<int32_t> types(nv, rD.mVarType);
        int32_t zone = 0;
        if (rZ.mType == 0)
            check(tecZoneCreateIJK(h, rZ.mTitle.c_str(), rZ.mI, rZ.mJ, rZ.mK, types.data(),
                                   rZ.mShare.data(), rZ.mLocation.data(), rZ.mPassive.data(),
                                   rZ.mConnShare, 0, 0, &zone),
                  "tecZoneCreateIJK");
        else
            check(tecZoneCreateFE(h, rZ.mTitle.c_str(), rZ.mType, rZ.mI, rZ.mJ, types.data(),
                                  rZ.mShare.data(), rZ.mLocation.data(), rZ.mPassive.data(),
                                  rZ.mConnShare, 0, 0, &zone),
                  "tecZoneCreateFE");
        if (rZ.mUnsteady)
            check(tecZoneSetUnsteadyOptions(h, zone, rZ.mTime, rZ.mStrand),
                  "tecZoneSetUnsteadyOptions");
        for (std::size_t v = 0; v < nv; ++v) {
            const std::vector<double>& rV = rZ.mValues[v];
            if (rV.empty())
                continue;
            const auto var = static_cast<int32_t>(v + 1);
            const auto n = static_cast<int64_t>(rV.size());
            if (rD.mVarType == 1) {
                const std::vector<float> f(rV.begin(), rV.end());
                check(tecZoneVarWriteFloatValues(h, zone, var, 0, n, f.data()), "float values");
            } else {
                check(tecZoneVarWriteDoubleValues(h, zone, var, 0, n, rV.data()), "double values");
            }
        }
        if (!rZ.mNodes.empty())
            check(tecZoneNodeMapWrite32(h, zone, 0, 1, static_cast<int64_t>(rZ.mNodes.size()),
                                        rZ.mNodes.data()),
                  "tecZoneNodeMapWrite32");
    }
    check(tecFileWriterClose(&h), "tecFileWriterClose");
}

void write_dat(const std::string& rName, const DataSet& rD) {
    std::ofstream out(rName);
    out.precision(17);
    out << "TITLE = \"meshio++ fixture\"\nVARIABLES =";
    for (const std::string& rV : rD.mVariables)
        out << " \"" << rV << "\"";
    out << "\n";
    const std::size_t nv = rD.mVariables.size();
    for (const Zone& rZ : rD.mZones) {
        out << "ZONE T=\"" << rZ.mTitle << "\", ZONETYPE=" << kZoneTypes[rZ.mType]
            << ", DATAPACKING=BLOCK";
        if (rZ.mType == 0)
            out << ", I=" << rZ.mI << ", J=" << rZ.mJ << ", K=" << rZ.mK;
        else
            out << ", NODES=" << rZ.mI << ", ELEMENTS=" << rZ.mJ;
        std::ostringstream cc, share, passive;
        for (std::size_t v = 0; v < nv; ++v) {
            if (rZ.mLocation[v] == 0)
                cc << (cc.tellp() > 0 ? "," : "") << v + 1;
            if (rZ.mShare[v] > 0)
                share << (share.tellp() > 0 ? "," : "") << "[" << v + 1 << "]=" << rZ.mShare[v];
            if (rZ.mPassive[v])
                passive << (passive.tellp() > 0 ? "," : "") << "[" << v + 1 << "]";
        }
        if (cc.tellp() > 0)
            out << ", VARLOCATION=([" << cc.str() << "]=CELLCENTERED)";
        if (share.tellp() > 0)
            out << ", VARSHARELIST=(" << share.str() << ")";
        if (passive.tellp() > 0)
            out << ", PASSIVEVARLIST=(" << passive.str() << ")";
        if (rZ.mConnShare > 0)
            out << ", CONNECTIVITYSHAREZONE=" << rZ.mConnShare;
        if (rZ.mUnsteady)
            out << ", SOLUTIONTIME=" << rZ.mTime << ", STRANDID=" << rZ.mStrand;
        out << "\n";
        for (std::size_t v = 0; v < nv; ++v) {
            if (rZ.mValues[v].empty())
                continue;
            for (std::size_t i = 0; i < rZ.mValues[v].size(); ++i)
                out << rZ.mValues[v][i] << ((i + 1) % 10 == 0 ? "\n" : " ");
            out << "\n";
        }
        if (!rZ.mNodes.empty()) {
            const std::size_t npc = rZ.mNodes.size() / static_cast<std::size_t>(rZ.mJ);
            for (std::size_t i = 0; i < rZ.mNodes.size(); ++i)
                out << rZ.mNodes[i] << ((i + 1) % npc == 0 ? "\n" : " ");
        }
    }
}

Zone zone(const std::string& rTitle, int32_t Type, int64_t I, int64_t J, int64_t K,
          std::vector<int32_t> Location) {
    Zone z;
    z.mTitle = rTitle;
    z.mType = Type;
    z.mI = I;
    z.mJ = J;
    z.mK = K;
    const std::size_t nv = Location.size();
    z.mLocation = std::move(Location);
    z.mShare.assign(nv, 0);
    z.mPassive.assign(nv, 0);
    z.mValues.assign(nv, {});
    return z;
}

// Two tetrahedra and two bricks in two FE zones; P nodal, Q cell-centred.
DataSet fe_mixed() {
    DataSet d;
    d.mVariables = {"X", "Y", "Z", "P", "Q"};
    Zone tets = zone("tets", 4, 5, 2, 0, {1, 1, 1, 1, 0});
    tets.mValues = {
        {0, 1, 0, 0, 1}, {0, 0, 1, 0, 1}, {0, 0, 0, 1, 1}, {0.5, 1.5, 2.5, 3.5, 4.5}, {10, 20}};
    tets.mNodes = {1, 2, 3, 4, 2, 5, 3, 4};
    Zone bricks = zone("bricks", 5, 12, 2, 0, {1, 1, 1, 1, 0});
    bricks.mValues.assign(5, {});
    for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 3; ++i) {
                bricks.mValues[0].push_back(2.0 + i);
                bricks.mValues[1].push_back(j);
                bricks.mValues[2].push_back(k);
                bricks.mValues[3].push_back(i + 10.0 * j + 100.0 * k);
            }
    bricks.mValues[4] = {-1, -2};
    // node (i, j, k) is 1 + i + 3 j + 6 k
    bricks.mNodes = {1, 2, 5, 4, 7, 8, 11, 10, 2, 3, 6, 5, 8, 9, 12, 11};
    d.mZones = {tets, bricks};
    return d;
}

// An ordered 3 x 2 x 2 zone; P nodal, Q cell-centred.
DataSet ordered() {
    DataSet d;
    d.mVariables = {"X", "Y", "Z", "P", "Q"};
    Zone z = zone("block", 0, 3, 2, 2, {1, 1, 1, 1, 0});
    for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 3; ++i) {
                z.mValues[0].push_back(0.5 * i);
                z.mValues[1].push_back(j);
                z.mValues[2].push_back(2.0 * k);
                z.mValues[3].push_back(i * j + k);
            }
    z.mValues[4] = {7, 8};
    d.mZones = {z};
    return d;
}

// Three steps of a quad zone (strand 1): the later zones share X, Y and the
// connectivity of the first; R is passive after the first step.
DataSet transient() {
    DataSet d;
    d.mVariables = {"X", "Y", "P", "R"};
    const double times[] = {0.0, 0.5, 1.0};
    for (int s = 0; s < 3; ++s) {
        Zone z = zone("step " + std::to_string(s), 3, 6, 2, 0, {1, 1, 1, 1});
        z.mUnsteady = true;
        z.mTime = times[s];
        z.mStrand = 1;
        for (int n = 0; n < 6; ++n)
            z.mValues[2].push_back(n + 10.0 * s);
        if (s == 0) {
            z.mValues[0] = {0, 1, 2, 0, 1, 2};
            z.mValues[1] = {0, 0, 0, 1, 1, 1};
            z.mValues[3] = {1, 1, 1, 2, 2, 2};
            z.mNodes = {1, 2, 5, 4, 2, 3, 6, 5};
        } else {
            z.mShare[0] = z.mShare[1] = 1;
            z.mPassive[3] = 1;
            z.mConnShare = 1;
        }
        d.mZones.push_back(z);
    }
    return d;
}

// Two-dimensional triangles with single-precision data (exact in float).
DataSet triangles() {
    DataSet d;
    d.mVariables = {"X", "Y", "T"};
    d.mVarType = 1;
    Zone z = zone("tris", 2, 4, 2, 0, {1, 1, 0});
    z.mValues = {{0, 1, 1, 0}, {0, 0, 1, 1}, {1.25, 2.5}};
    z.mNodes = {1, 2, 3, 1, 3, 4};
    d.mZones = {z};
    return d;
}

}  // namespace

int main() {
    const std::pair<const char*, DataSet> sets[] = {
        {"fe_mixed", fe_mixed()},
        {"ordered", ordered()},
        {"transient", transient()},
        {"triangles", triangles()},
    };
    for (const auto& [name, data] : sets) {
        write_szplt(std::string(name) + ".szplt", data);
        write_dat(std::string(name) + ".dat", data);
    }
    return 0;
}
