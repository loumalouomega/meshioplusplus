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
// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/dxf.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "../detail/open_source.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

namespace {

using DxfV3 = std::array<double, 3>;
constexpr double kDxfPi = 3.14159265358979323846;
constexpr int kDxfCircleSegments = 64;
constexpr std::array<std::pair<int, double>, 11> kDxfUnits{{{1, 0.0254},
                                                            {2, 0.3048},
                                                            {3, 1609.344},
                                                            {4, 0.001},
                                                            {5, 0.01},
                                                            {6, 1.0},
                                                            {7, 1000.0},
                                                            {10, 0.9144},
                                                            {12, 1e-9},
                                                            {13, 1e-6},
                                                            {14, 0.1}}};

DxfV3 dxf_cross(const DxfV3& rA, const DxfV3& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}
double dxf_dot(const DxfV3& rA, const DxfV3& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}
DxfV3 dxf_unit(const DxfV3& rA) {
    const double n = std::sqrt(rA[0] * rA[0] + rA[1] * rA[1] + rA[2] * rA[2]);
    return {rA[0] / n, rA[1] / n, rA[2] / n};
}

// The DXF arbitrary-axis algorithm: the object coordinate system of `rN`.
std::array<DxfV3, 3> dxf_ocs(const DxfV3& rN) {
    const DxfV3 n = dxf_unit(rN);
    DxfV3 ax = (std::fabs(n[0]) < 1.0 / 64.0 && std::fabs(n[1]) < 1.0 / 64.0)
                   ? dxf_cross({0.0, 1.0, 0.0}, n)
                   : dxf_cross({0.0, 0.0, 1.0}, n);
    ax = dxf_unit(ax);
    return {ax, dxf_cross(n, ax), n};
}

std::string dxf_num(double V) {
    if (V == 0.0)
        V = 0.0;  // no "-0"
    char buf[40];
    detail::snprintf_c(buf, sizeof(buf), "%.17g", V);
    return buf;
}

class DxfOut {
public:
    explicit DxfOut(int Handle) : mHandle(Handle) {}
    void Pair(int Code, const std::string& rValue) {
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%3d", Code);
        mText += buf;
        mText += '\n';
        mText += rValue;
        mText += '\n';
    }
    void Pair(int Code, long long Value) { Pair(Code, std::to_string(Value)); }
    std::string Next() {
        char buf[24];
        std::snprintf(buf, sizeof(buf), "%X", static_cast<unsigned>(mHandle));
        ++mHandle;
        return buf;
    }
    int Handle() const { return mHandle; }
    std::string mText;

private:
    int mHandle;
};

int dxf_insunits(const Mesh& rMesh) {
    if (!rMesh.HasFieldData(detail::kLengthUnitName))
        return 0;
    const NDArray& v = rMesh.FieldData(detail::kLengthUnitName);
    if (v.Size() != 1)
        return 0;
    const double m = detail::read_double(v, 0);
    for (const auto& [code, metres] : kDxfUnits)
        if (std::fabs(m - metres) <= 1e-9 * metres)
            return code;
    return 0;
}

void dxf_entity_head(DxfOut& rOut, const char* pType, const char* pSub) {
    rOut.Pair(0, pType);
    rOut.Pair(5, rOut.Next());
    rOut.Pair(330, "1F");
    rOut.Pair(100, "AcDbEntity");
    rOut.Pair(8, "0");
    rOut.Pair(100, pSub);
}

void dxf_entities(DxfOut& rOut, const Mesh& rMesh) {
    const bool sketch = rMesh.HasFieldData("sketch:type") && rMesh.HasFieldData("sketch:params") &&
                        rMesh.HasFieldData("sketch:frame") && rMesh.HasFieldData("sketch:frames");
    if (sketch) {
        const NDArray& types = rMesh.FieldData("sketch:type");
        const NDArray& params = rMesh.FieldData("sketch:params");
        const NDArray& frame_of = rMesh.FieldData("sketch:frame");
        const NDArray& frames = rMesh.FieldData("sketch:frames");
        if (params.Size() != 10 * types.Size() || frame_of.Size() != types.Size())
            throw WriteError("meshio++: dxf: the sketch:* field data do not match");
        for (std::size_t k = 0; k < types.Size(); ++k) {
            const std::int64_t t = detail::read_int(types, k);
            const std::int64_t f = detail::read_int(frame_of, k);
            if (f < 0 || static_cast<std::size_t>(16 * (f + 1)) > frames.Size())
                throw WriteError("meshio++: dxf: a sketch piece names a missing frame");
            const auto m = [&](std::size_t R, std::size_t C) {
                return detail::read_double(frames, static_cast<std::size_t>(16 * f) + 4 * R + C);
            };
            DxfV3 e1, e2, n, o;
            for (std::size_t r = 0; r < 3; ++r) {
                e1[r] = m(r, 0);
                e2[r] = m(r, 1);
                n[r] = m(r, 2);
                o[r] = m(r, 3);
            }
            double p[10];
            for (std::size_t q = 0; q < 10; ++q)
                p[q] = detail::read_double(params, 10 * k + q);
            const auto world = [&](double U, double V) {
                DxfV3 w;
                for (std::size_t r = 0; r < 3; ++r)
                    w[r] = o[r] + U * e1[r] + V * e2[r];
                return w;
            };
            const DxfV3 start = world(p[0], p[1]), end = world(p[2], p[3]),
                        centre = world(p[4], p[5]);
            if (t == 1) {
                dxf_entity_head(rOut, "LINE", "AcDbLine");
                for (std::size_t r = 0; r < 3; ++r)
                    rOut.Pair(10 + 10 * static_cast<int>(r), dxf_num(start[r]));
                for (std::size_t r = 0; r < 3; ++r)
                    rOut.Pair(11 + 10 * static_cast<int>(r), dxf_num(end[r]));
                continue;
            }
            const std::array<DxfV3, 3> ocs = dxf_ocs(n);
            const DxfV3 c{dxf_dot(centre, ocs[0]), dxf_dot(centre, ocs[1]),
                          dxf_dot(centre, ocs[2])};
            dxf_entity_head(rOut, t == 3 ? "CIRCLE" : "ARC", "AcDbCircle");
            for (std::size_t r = 0; r < 3; ++r)
                rOut.Pair(10 + 10 * static_cast<int>(r), dxf_num(c[r]));
            rOut.Pair(40, dxf_num(p[6]));
            for (std::size_t r = 0; r < 3; ++r)
                rOut.Pair(210 + 10 * static_cast<int>(r), dxf_num(ocs[2][r]));
            if (t == 2) {
                const DxfV3 ds{start[0] - centre[0], start[1] - centre[1], start[2] - centre[2]};
                const DxfV3 de{end[0] - centre[0], end[1] - centre[1], end[2] - centre[2]};
                double a0 = std::atan2(dxf_dot(ds, ocs[1]), dxf_dot(ds, ocs[0])) * 180.0 / kDxfPi;
                double a1 = std::atan2(dxf_dot(de, ocs[1]), dxf_dot(de, ocs[0])) * 180.0 / kDxfPi;
                if (p[7] < 0.0)
                    std::swap(a0, a1);
                rOut.Pair(100, "AcDbArc");
                rOut.Pair(50, dxf_num(a0));
                rOut.Pair(51, dxf_num(a1));
            }
        }
        return;
    }
    const std::size_t dim = rMesh.PointDim();
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cells = rMesh.Cells(b);
        if (cells.Type() == "vertex")
            continue;
        if (cells.Type() != "line")
            throw WriteError("meshio++: dxf: only line cells or a sketch are written, not " +
                             cells.Type() + " (slice or sections give lines; fit_curves a sketch)");
        const NDArray& conn = cells.Conn();
        for (std::size_t s = 0; s < cells.NumCells(); ++s) {
            dxf_entity_head(rOut, "LINE", "AcDbLine");
            for (int e = 0; e < 2; ++e) {
                const auto v = static_cast<std::size_t>(detail::read_int(conn, 2 * s + e));
                for (std::size_t r = 0; r < 3; ++r)
                    rOut.Pair(
                        (e == 0 ? 10 : 11) + 10 * static_cast<int>(r),
                        dxf_num(r < dim ? detail::read_double(rMesh.Points(), v * dim + r) : 0.0));
            }
        }
    }
}

void dxf_table(DxfOut& rOut, const char* pName, const char* pHandle, int Count) {
    rOut.Pair(0, "TABLE");
    rOut.Pair(2, pName);
    rOut.Pair(5, pHandle);
    rOut.Pair(330, "0");
    rOut.Pair(100, "AcDbSymbolTable");
    rOut.Pair(70, Count);
}

void dxf_record(DxfOut& rOut, const char* pType, const char* pHandle, const char* pOwner,
                const char* pSub) {
    rOut.Pair(0, pType);
    rOut.Pair(5, pHandle);
    rOut.Pair(330, pOwner);
    rOut.Pair(100, "AcDbSymbolTableRecord");
    rOut.Pair(100, pSub);
}

std::string dxf_document(const Mesh& rMesh) {
    DxfOut body(0x30);
    dxf_entities(body, rMesh);
    DxfOut out(0x10);
    out.Pair(0, "SECTION");
    out.Pair(2, "HEADER");
    out.Pair(9, "$ACADVER");
    out.Pair(1, "AC1015");
    out.Pair(9, "$HANDSEED");
    {
        char buf[24];
        std::snprintf(buf, sizeof(buf), "%X", static_cast<unsigned>(body.Handle()));
        out.Pair(5, buf);
    }
    out.Pair(9, "$INSUNITS");
    out.Pair(70, dxf_insunits(rMesh));
    out.Pair(0, "ENDSEC");
    out.Pair(0, "SECTION");
    out.Pair(2, "TABLES");
    dxf_table(out, "VPORT", "8", 0);
    out.Pair(0, "ENDTAB");
    dxf_table(out, "LTYPE", "5", 3);
    const std::array<std::array<const char*, 3>, 3> ltypes{
        {{"ByBlock", "14", ""}, {"ByLayer", "15", ""}, {"Continuous", "16", "Solid line"}}};
    for (const auto& r_l : ltypes) {
        dxf_record(out, "LTYPE", r_l[1], "5", "AcDbLinetypeTableRecord");
        out.Pair(2, r_l[0]);
        out.Pair(70, 0);
        out.Pair(3, r_l[2]);
        out.Pair(72, 65);
        out.Pair(73, 0);
        out.Pair(40, "0.0");
    }
    out.Pair(0, "ENDTAB");
    dxf_table(out, "LAYER", "2", 1);
    dxf_record(out, "LAYER", "10", "2", "AcDbLayerTableRecord");
    out.Pair(2, "0");
    out.Pair(70, 0);
    out.Pair(62, 7);
    out.Pair(6, "Continuous");
    out.Pair(0, "ENDTAB");
    dxf_table(out, "STYLE", "3", 1);
    dxf_record(out, "STYLE", "11", "3", "AcDbTextStyleTableRecord");
    out.Pair(2, "Standard");
    out.Pair(70, 0);
    out.Pair(40, "0.0");
    out.Pair(41, "1.0");
    out.Pair(50, "0.0");
    out.Pair(71, 0);
    out.Pair(42, "2.5");
    out.Pair(3, "txt");
    out.Pair(4, "");
    out.Pair(0, "ENDTAB");
    dxf_table(out, "VIEW", "6", 0);
    out.Pair(0, "ENDTAB");
    dxf_table(out, "UCS", "7", 0);
    out.Pair(0, "ENDTAB");
    dxf_table(out, "APPID", "9", 1);
    dxf_record(out, "APPID", "12", "9", "AcDbRegAppTableRecord");
    out.Pair(2, "ACAD");
    out.Pair(70, 0);
    out.Pair(0, "ENDTAB");
    dxf_table(out, "DIMSTYLE", "A", 1);
    out.Pair(100, "AcDbDimStyleTable");
    out.Pair(0, "DIMSTYLE");
    out.Pair(105, "13");
    out.Pair(330, "A");
    out.Pair(100, "AcDbSymbolTableRecord");
    out.Pair(100, "AcDbDimStyleTableRecord");
    out.Pair(2, "Standard");
    out.Pair(70, 0);
    out.Pair(0, "ENDTAB");
    dxf_table(out, "BLOCK_RECORD", "1", 2);
    for (const auto& r_b : std::array<std::array<const char*, 3>, 2>{
             {{"*Model_Space", "1F", "22"}, {"*Paper_Space", "1B", "1E"}}}) {
        dxf_record(out, "BLOCK_RECORD", r_b[1], "1", "AcDbBlockTableRecord");
        out.Pair(2, r_b[0]);
        out.Pair(340, r_b[2]);
    }
    out.Pair(0, "ENDTAB");
    out.Pair(0, "ENDSEC");
    out.Pair(0, "SECTION");
    out.Pair(2, "BLOCKS");
    for (const auto& r_b : std::array<std::array<const char*, 4>, 2>{
             {{"*Model_Space", "1F", "20", "21"}, {"*Paper_Space", "1B", "1C", "1D"}}}) {
        out.Pair(0, "BLOCK");
        out.Pair(5, r_b[2]);
        out.Pair(330, r_b[1]);
        out.Pair(100, "AcDbEntity");
        out.Pair(8, "0");
        out.Pair(100, "AcDbBlockBegin");
        out.Pair(2, r_b[0]);
        out.Pair(70, 0);
        out.Pair(10, "0.0");
        out.Pair(20, "0.0");
        out.Pair(30, "0.0");
        out.Pair(3, r_b[0]);
        out.Pair(1, "");
        out.Pair(0, "ENDBLK");
        out.Pair(5, r_b[3]);
        out.Pair(330, r_b[1]);
        out.Pair(100, "AcDbEntity");
        out.Pair(8, "0");
        out.Pair(100, "AcDbBlockEnd");
    }
    out.Pair(0, "ENDSEC");
    out.Pair(0, "SECTION");
    out.Pair(2, "ENTITIES");
    out.mText += body.mText;
    out.Pair(0, "ENDSEC");
    out.Pair(0, "SECTION");
    out.Pair(2, "OBJECTS");
    struct DxfDict {
        const char* mHandle;
        const char* mOwner;
        std::vector<std::pair<const char*, const char*>> mEntries;
    };
    const std::array<DxfDict, 3> dicts{{{"C", "0", {{"ACAD_GROUP", "D"}, {"ACAD_LAYOUT", "1A"}}},
                                        {"D", "C", {}},
                                        {"1A", "C", {{"Layout1", "1E"}, {"Model", "22"}}}}};
    for (const DxfDict& r_d : dicts) {
        out.Pair(0, "DICTIONARY");
        out.Pair(5, r_d.mHandle);
        out.Pair(330, r_d.mOwner);
        out.Pair(100, "AcDbDictionary");
        out.Pair(281, 1);
        for (const auto& r_e : r_d.mEntries) {
            out.Pair(3, r_e.first);
            out.Pair(350, r_e.second);
        }
    }
    for (const auto& r_l : std::array<std::array<const char*, 4>, 2>{
             {{"22", "Model", "1F", "0"}, {"1E", "Layout1", "1B", "1"}}}) {
        out.Pair(0, "LAYOUT");
        out.Pair(5, r_l[0]);
        out.Pair(330, "1A");
        out.Pair(100, "AcDbPlotSettings");
        out.Pair(1, "");
        out.Pair(4, "");
        out.Pair(6, "");
        for (const int code : {40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 140, 141})
            out.Pair(code, "0.0");
        out.Pair(142, "1.0");
        out.Pair(143, "1.0");
        out.Pair(70, 688);
        out.Pair(72, 0);
        out.Pair(73, 0);
        out.Pair(74, 5);
        out.Pair(7, "");
        out.Pair(75, 16);
        out.Pair(147, "1.0");
        out.Pair(148, "0.0");
        out.Pair(149, "0.0");
        out.Pair(100, "AcDbLayout");
        out.Pair(1, r_l[1]);
        out.Pair(70, 1);
        out.Pair(71, r_l[3]);
        out.Pair(10, "0.0");
        out.Pair(20, "0.0");
        out.Pair(11, "420.0");
        out.Pair(21, "297.0");
        for (const int code : {12, 22, 32, 14, 24, 34, 15, 25, 35, 146, 13, 23, 33})
            out.Pair(code, "0.0");
        out.Pair(16, "1.0");
        out.Pair(26, "0.0");
        out.Pair(36, "0.0");
        out.Pair(17, "0.0");
        out.Pair(27, "1.0");
        out.Pair(37, "0.0");
        out.Pair(76, 0);
        out.Pair(330, r_l[2]);
    }
    out.Pair(0, "ENDSEC");
    out.Pair(0, "EOF");
    return out.mText;
}

// ---- reading -------------------------------------------------------------------

std::string dxf_strip(const std::string& rS) {
    const std::size_t a = rS.find_first_not_of(" \t\r\n\v\f");
    if (a == std::string::npos)
        return {};
    const std::size_t b = rS.find_last_not_of(" \t\r\n\v\f");
    return rS.substr(a, b - a + 1);
}

double dxf_float(const std::string& rValue, int Code) {
    const char* begin = rValue.c_str();
    const char* end = begin;
    const double v = detail::parse_double(begin, end);
    if (end == begin || *end != '\0')
        throw ReadError("meshio++: dxf: group " + std::to_string(Code) + ": '" + rValue +
                        "' is not a number");
    if (!std::isfinite(v))
        throw ReadError("meshio++: dxf: group " + std::to_string(Code) + ": not finite");
    return v;
}

struct DxfReader {
    std::vector<DxfV3> mPoints;
    std::map<DxfV3, std::int64_t> mIndex;
    std::vector<std::int64_t> mConn, mEntity, mType;

    std::int64_t Point(const DxfV3& rP) {
        const auto it = mIndex.find(rP);
        if (it != mIndex.end())
            return it->second;
        const auto id = static_cast<std::int64_t>(mPoints.size());
        mIndex.emplace(rP, id);
        mPoints.push_back(rP);
        return id;
    }
    void Polyline(const std::vector<DxfV3>& rPts, bool Closed, std::int64_t Entity, int Kind) {
        std::vector<std::int64_t> ids;
        for (const DxfV3& r_p : rPts)
            ids.push_back(Point(r_p));
        if (Closed && !ids.empty())
            ids.push_back(ids.front());
        for (std::size_t k = 0; k + 1 < ids.size(); ++k)
            if (ids[k] != ids[k + 1]) {
                mConn.push_back(ids[k]);
                mConn.push_back(ids[k + 1]);
                mEntity.push_back(Entity);
                mType.push_back(Kind);
            }
    }
};

std::vector<DxfV3> dxf_arc_points(const DxfV3& rC, double R, double A0, double A1,
                                  const std::array<DxfV3, 3>& rOcs, bool Closed) {
    double sweep = 360.0;
    if (!Closed) {
        sweep = std::fmod(A1 - A0, 360.0);
        if (sweep < 0.0)
            sweep += 360.0;
        if (sweep == 0.0)
            sweep = 360.0;
    }
    const int k = std::max(4, static_cast<int>(std::ceil(sweep / 360.0 * kDxfCircleSegments)));
    const int count = Closed ? k : k + 1;
    std::vector<DxfV3> out;
    for (int i = 0; i < count; ++i) {
        const double a = (A0 + sweep * i / k) * kDxfPi / 180.0;
        const double lx = rC[0] + R * std::cos(a), ly = rC[1] + R * std::sin(a);
        DxfV3 w;
        for (std::size_t r = 0; r < 3; ++r)
            w[r] = lx * rOcs[0][r] + ly * rOcs[1][r] + rC[2] * rOcs[2][r];
        out.push_back(w);
    }
    return out;
}

std::vector<DxfV3> dxf_bspline(int P, const std::vector<double>& rKnots,
                               const std::vector<DxfV3>& rCtrl, const std::vector<double>& rW,
                               int N) {
    std::vector<std::array<double, 4>> hom(rCtrl.size());
    for (std::size_t i = 0; i < rCtrl.size(); ++i) {
        const double w = rW.empty() ? 1.0 : rW[i];
        hom[i] = {rCtrl[i][0] * w, rCtrl[i][1] * w, rCtrl[i][2] * w, w};
    }
    const double t0 = rKnots[static_cast<std::size_t>(P)];
    const double t1 = rKnots[rKnots.size() - static_cast<std::size_t>(P) - 1];
    std::vector<DxfV3> out;
    for (int i = 0; i <= N; ++i) {
        const double t = t0 + (t1 - t0) * i / N;
        std::size_t k = static_cast<std::size_t>(P);
        while (k < rKnots.size() - static_cast<std::size_t>(P) - 2 && rKnots[k + 1] <= t)
            ++k;
        std::vector<std::array<double, 4>> d(static_cast<std::size_t>(P) + 1);
        for (std::size_t j = 0; j <= static_cast<std::size_t>(P); ++j)
            d[j] = hom[j + k - static_cast<std::size_t>(P)];
        for (std::size_t r = 1; r <= static_cast<std::size_t>(P); ++r)
            for (std::size_t j = static_cast<std::size_t>(P); j >= r; --j) {
                const double den =
                    rKnots[j + 1 + k - r] - rKnots[j + k - static_cast<std::size_t>(P)];
                const double alpha =
                    den == 0.0 ? 0.0 : (t - rKnots[j + k - static_cast<std::size_t>(P)]) / den;
                for (std::size_t c = 0; c < 4; ++c)
                    d[j][c] = (1.0 - alpha) * d[j - 1][c] + alpha * d[j][c];
            }
        const auto& h = d[static_cast<std::size_t>(P)];
        out.push_back({h[0] / h[3], h[1] / h[3], h[2] / h[3]});
    }
    return out;
}

}  // namespace

void write_dxf(const std::string& rPath, const Mesh& rMesh) {
    // The provenance block leads as 999 comment pairs, which every DXF reader skips.
    const std::string text =
        detail::provenance_render_lines(detail::SlotTier::Block, "999\n") + dxf_document(rMesh);
    auto out = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!out)
        throw WriteError("meshio++: dxf: cannot open '" + rPath + "' for writing");
    out << text;
    if (!out)
        throw WriteError("meshio++: dxf: cannot write '" + rPath + "'");
}

Mesh read_dxf(const std::string& rPath) {
    const detail::FileSource source =
        detail::open_source(rPath, "meshio++: dxf: cannot open '" + rPath + "'");
    const std::string_view text(source.Data(), source.Size());
    std::vector<std::string> lines;
    {
        std::size_t a = 0;
        while (a < text.size()) {
            std::size_t b = text.find('\n', a);
            if (b == std::string::npos)
                b = text.size();
            lines.emplace_back(text.substr(a, b - a));
            a = b + 1;
        }
    }
    if (lines.size() % 2 == 1 && dxf_strip(lines.back()).empty())
        lines.pop_back();
    if (lines.size() % 2 == 1)
        throw ReadError("meshio++: dxf: the file ends inside a group code pair");
    std::vector<std::pair<int, std::string>> pairs;
    pairs.reserve(lines.size() / 2);
    for (std::size_t k = 0; k < lines.size(); k += 2) {
        const std::string code = dxf_strip(lines[k]);
        const char* begin = code.c_str();
        const char* end = begin;
        const double c = detail::parse_double(begin, end);
        if (end == begin || *end != '\0' || c != std::floor(c) || std::fabs(c) > 1e6)
            throw ReadError("meshio++: dxf: line " + std::to_string(k + 1) + ": '" + code +
                            "' is not a group code");
        pairs.emplace_back(static_cast<int>(c), dxf_strip(lines[k + 1]));
    }
    DxfReader rd;
    int units = 0;
    std::string section;
    std::int64_t entity = 0;
    std::size_t k = 0;
    while (k < pairs.size()) {
        const auto& [code, value] = pairs[k];
        if (code == 0 && value == "SECTION" && k + 1 < pairs.size()) {
            section = pairs[k + 1].second;
            k += 2;
            continue;
        }
        if (code == 0 && value == "ENDSEC") {
            section.clear();
            ++k;
            continue;
        }
        if (section == "HEADER" && code == 9 && value == "$INSUNITS" && k + 1 < pairs.size()) {
            units = static_cast<int>(dxf_float(pairs[k + 1].second, 70));
            k += 2;
            continue;
        }
        if (section != "ENTITIES" || code != 0) {
            ++k;
            continue;
        }
        const std::string kind = value;
        if (kind != "LINE" && kind != "ARC" && kind != "CIRCLE" && kind != "LWPOLYLINE" &&
            kind != "SPLINE")
            throw ReadError("meshio++: dxf: the entity " + kind +
                            " is not read (only LINE, ARC, CIRCLE, LWPOLYLINE and SPLINE)");
        ++k;
        std::vector<std::pair<int, std::string>> groups;
        while (k < pairs.size() && pairs[k].first != 0)
            groups.push_back(pairs[k++]);
        std::map<int, std::vector<std::string>> g;
        for (const auto& [c, v] : groups)
            g[c].push_back(v);
        const auto num = [&](int C, double Default) {
            const auto it = g.find(C);
            return it == g.end() ? Default : dxf_float(it->second.front(), C);
        };
        const DxfV3 n{num(210, 0.0), num(220, 0.0), num(230, 1.0)};
        if (n[0] == 0.0 && n[1] == 0.0 && n[2] == 0.0)
            throw ReadError("meshio++: dxf: " + kind + ": zero extrusion direction");
        const std::array<DxfV3, 3> ocs = dxf_ocs(n);
        if (kind == "LINE") {
            rd.Polyline({{num(10, 0.0), num(20, 0.0), num(30, 0.0)},
                         {num(11, 0.0), num(21, 0.0), num(31, 0.0)}},
                        false, entity, 1);
        } else if (kind == "ARC" || kind == "CIRCLE") {
            const double r = num(40, 0.0);
            if (!(r > 0.0))
                throw ReadError("meshio++: dxf: " + kind + ": radius must be positive");
            const DxfV3 c{num(10, 0.0), num(20, 0.0), num(30, 0.0)};
            const bool closed = kind == "CIRCLE";
            const double a0 = closed ? 0.0 : num(50, 0.0), a1 = closed ? 360.0 : num(51, 0.0);
            rd.Polyline(dxf_arc_points(c, r, a0, a1, ocs, closed), closed, entity,
                        kind == "ARC" ? 2 : 3);
        } else if (kind == "LWPOLYLINE") {
            const std::vector<std::string> xs = g.count(10) ? g[10] : std::vector<std::string>{};
            const std::vector<std::string> ys = g.count(20) ? g[20] : std::vector<std::string>{};
            if (xs.size() != ys.size())
                throw ReadError("meshio++: dxf: LWPOLYLINE: unpaired vertex coordinates");
            std::vector<double> bulges(xs.size(), 0.0);
            std::size_t seen = 0;
            for (const auto& [c, v] : groups) {
                if (c == 10)
                    ++seen;
                else if (c == 42 && seen > 0)
                    bulges[seen - 1] = dxf_float(v, 42);
            }
            const double elev = num(38, 0.0);
            const bool closed = (static_cast<long long>(num(70, 0.0)) & 1) == 1;
            std::vector<std::array<double, 2>> verts;
            for (std::size_t i = 0; i < xs.size(); ++i)
                verts.push_back({dxf_float(xs[i], 10), dxf_float(ys[i], 20)});
            std::vector<std::array<double, 2>> local;
            const std::size_t count = verts.size();
            const std::size_t spans = closed ? count : (count > 0 ? count - 1 : 0);
            for (std::size_t i = 0; i < spans; ++i) {
                const auto& p0 = verts[i];
                const auto& p1 = verts[(i + 1) % count];
                local.push_back(p0);
                const double b = bulges[i];
                if (b != 0.0) {
                    const double theta = 4.0 * std::atan(b);
                    const double chord = std::hypot(p1[0] - p0[0], p1[1] - p0[1]);
                    if (chord > 0.0) {
                        const double r = chord / (2.0 * std::sin(std::fabs(theta) / 2.0));
                        const double mx = (p0[0] + p1[0]) / 2.0, my = (p0[1] + p1[1]) / 2.0;
                        const double h = r * std::cos(theta / 2.0) * (b > 0 ? 1 : -1);
                        const double dx = (p1[0] - p0[0]) / chord, dy = (p1[1] - p0[1]) / chord;
                        const double cx = mx - h * dy, cy = my + h * dx;
                        const double a0 = std::atan2(p0[1] - cy, p0[0] - cx);
                        const int steps =
                            std::max(2, static_cast<int>(std::ceil(std::fabs(theta) / (2 * kDxfPi) *
                                                                   kDxfCircleSegments)));
                        for (int s = 1; s < steps; ++s) {
                            const double a = a0 + theta * s / steps;
                            local.push_back({cx + r * std::cos(a), cy + r * std::sin(a)});
                        }
                    }
                }
            }
            if (!closed && count > 0)
                local.push_back(verts.back());
            std::vector<DxfV3> world;
            for (const auto& r_l : local) {
                DxfV3 w;
                for (std::size_t r = 0; r < 3; ++r)
                    w[r] = r_l[0] * ocs[0][r] + r_l[1] * ocs[1][r] + elev * ocs[2][r];
                world.push_back(w);
            }
            rd.Polyline(world, closed, entity, 4);
        } else {
            const int degree = static_cast<int>(num(71, 3.0));
            std::vector<double> knots, weights;
            for (const std::string& v : g.count(40) ? g[40] : std::vector<std::string>{})
                knots.push_back(dxf_float(v, 40));
            for (const std::string& v : g.count(41) ? g[41] : std::vector<std::string>{})
                weights.push_back(dxf_float(v, 41));
            const std::vector<std::string> xs = g.count(10) ? g[10] : std::vector<std::string>{};
            const std::vector<std::string> ys = g.count(20) ? g[20] : std::vector<std::string>{};
            const std::vector<std::string> zs = g.count(30) ? g[30] : std::vector<std::string>{};
            std::vector<DxfV3> ctrl;
            for (std::size_t i = 0; i < std::min({xs.size(), ys.size(), zs.size()}); ++i)
                ctrl.push_back({dxf_float(xs[i], 10), dxf_float(ys[i], 20), dxf_float(zs[i], 30)});
            if (degree < 1 || ctrl.size() <= static_cast<std::size_t>(degree) ||
                knots.size() != ctrl.size() + static_cast<std::size_t>(degree) + 1)
                throw ReadError("meshio++: dxf: SPLINE: inconsistent degree, knots and points");
            if (!weights.empty() && weights.size() != ctrl.size())
                throw ReadError("meshio++: dxf: SPLINE: one weight per control point");
            rd.Polyline(
                dxf_bspline(degree, knots, ctrl, weights, 8 * static_cast<int>(ctrl.size())), false,
                entity, 5);
        }
        ++entity;
    }
    Mesh mesh;
    NDArray pts = NDArray::Uninit(DType::Float64, {rd.mPoints.size(), std::size_t{3}});
    for (std::size_t i = 0; i < rd.mPoints.size(); ++i)
        for (std::size_t r = 0; r < 3; ++r)
            pts.As<double>()[3 * i + r] = rd.mPoints[i][r];
    mesh.AssignPoints(std::move(pts));
    if (!rd.mConn.empty()) {
        const std::size_t ns = rd.mEntity.size();
        NDArray c = NDArray::Uninit(DType::Int64, {ns, std::size_t{2}});
        std::copy(rd.mConn.begin(), rd.mConn.end(), c.As<std::int64_t>());
        mesh.AddCellBlock("line", std::move(c));
        NDArray e = NDArray::Uninit(DType::Int64, {ns});
        std::copy(rd.mEntity.begin(), rd.mEntity.end(), e.As<std::int64_t>());
        NDArray t = NDArray::Uninit(DType::Int64, {ns});
        std::copy(rd.mType.begin(), rd.mType.end(), t.As<std::int64_t>());
        std::vector<NDArray> eb, tb;
        eb.push_back(std::move(e));
        tb.push_back(std::move(t));
        mesh.AddCellData("dxf:entity", std::move(eb));
        mesh.AddCellData("dxf:type", std::move(tb));
    }
    for (const auto& [code, metres] : kDxfUnits)
        if (code == units) {
            NDArray u = NDArray::Uninit(DType::Float64, {std::size_t{1}});
            u.As<double>()[0] = metres;
            mesh.AddFieldData(detail::kLengthUnitName, std::move(u));
        }
    return mesh;
}

}  // namespace meshioplusplus
