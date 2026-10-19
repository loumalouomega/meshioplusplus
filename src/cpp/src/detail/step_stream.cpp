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

// The AP214 frame of a STEP B-rep. See detail/step_stream.hpp.

// System includes
#include <cstdint>
#include <ios>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"

// Project includes (private, not installed)
#include "step_stream.hpp"

namespace meshioplusplus {
namespace detail {
namespace {

std::string ss_basename(const std::string& rPath) {
    const std::size_t slash = rPath.find_last_of("/\\");
    return slash == std::string::npos ? rPath : rPath.substr(slash + 1);
}

}  // namespace

std::string step_real(double V) {
    char buf[64];
    snprintf_c(buf, sizeof(buf), "%.17g", V);
    std::string s(buf);
    const std::size_t e = s.find('e');
    std::string mant = e == std::string::npos ? s : s.substr(0, e);
    std::string expo = e == std::string::npos ? "" : "E" + s.substr(e + 1);
    if (mant.find('.') == std::string::npos)
        mant += '.';
    return mant + expo;
}

std::string step_string(const std::string& rText) {
    std::string out;
    for (unsigned char c : rText) {
        if (c == '\'')
            out += "''";
        else if (c >= 0x20 && c < 0x7f)
            out += static_cast<char>(c);
        else
            out += '_';
    }
    return out;
}

StepStream::StepStream(const std::string& rPath, const std::string& rDescription) : mPath(rPath) {
    const std::string name = step_string(ss_basename(rPath));
    mOut += "ISO-10303-21;\n/*\n";
    for (std::string line : provenance_lines(SlotTier::Block)) {
        for (std::size_t at = line.find("*/"); at != std::string::npos; at = line.find("*/", at))
            line.replace(at, 2, "* /");
        mOut += line + "\n";
    }
    mOut += "*/\nHEADER;\n";
    mOut += "FILE_DESCRIPTION(('" + step_string(rDescription) + "'),'2;1');\n";
    mOut += "FILE_NAME('" + name + "','1970-01-01T00:00:00',(''),(''),'meshio++','meshio++','');\n";
    mOut += "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n";

    const std::int64_t app =
        Ent("APPLICATION_CONTEXT('core data for automotive mechanical design processes')");
    Ent("APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000," +
        Ref(app) + ")");
    const std::int64_t pctx = Ent("PRODUCT_CONTEXT(''," + Ref(app) + ",'mechanical')");
    const std::int64_t pdctx =
        Ent("PRODUCT_DEFINITION_CONTEXT('part definition'," + Ref(app) + ",'design')");
    const std::int64_t prod =
        Ent("PRODUCT('" + name + "','" + name + "',''," + "(" + Ref(pctx) + "))");
    Ent("PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(" + Ref(prod) + "))");
    const std::int64_t form = Ent("PRODUCT_DEFINITION_FORMATION('',''," + Ref(prod) + ")");
    const std::int64_t pdef =
        Ent("PRODUCT_DEFINITION('design',''," + Ref(form) + "," + Ref(pdctx) + ")");
    mShape = Ent("PRODUCT_DEFINITION_SHAPE('',''," + Ref(pdef) + ")");
    const std::int64_t len_unit = Ent("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))");
    const std::int64_t ang_unit = Ent("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))");
    const std::int64_t sol_unit = Ent("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())");
    const std::int64_t unc =
        Ent("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07)," + Ref(len_unit) +
            ",'distance_accuracy_value','confusion accuracy')");
    mContext =
        Ent("(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((" + Ref(unc) +
            "))GLOBAL_UNIT_ASSIGNED_CONTEXT((" + Ref(len_unit) + "," + Ref(ang_unit) + "," +
            Ref(sol_unit) +
            "))REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY'))");
    const std::int64_t o_pt = Ent("CARTESIAN_POINT('',(0.,0.,0.))");
    const std::int64_t o_z = Ent("DIRECTION('',(0.,0.,1.))");
    const std::int64_t o_x = Ent("DIRECTION('',(1.,0.,0.))");
    mOrigin = Ent("AXIS2_PLACEMENT_3D(''," + Ref(o_pt) + "," + Ref(o_z) + "," + Ref(o_x) + ")");
}

std::int64_t StepStream::Ent(const std::string& rBody) {
    const std::int64_t id = mNext++;
    mOut += '#';
    mOut += std::to_string(id);
    mOut += '=';
    mOut += rBody;
    mOut += ";\n";
    return id;
}

std::string StepStream::Ref(std::int64_t Id) {
    return "#" + std::to_string(Id);
}

std::string StepStream::Triple(const double* pV) {
    return "(" + step_real(pV[0]) + "," + step_real(pV[1]) + "," + step_real(pV[2]) + ")";
}

void StepStream::Finish(const std::vector<std::int64_t>& rSolids) {
    std::string solids;
    for (std::int64_t s : rSolids)
        solids += "," + Ref(s);
    const std::int64_t rep = Ent("ADVANCED_BREP_SHAPE_REPRESENTATION('',(" + Ref(mOrigin) + solids +
                                 ")," + Ref(mContext) + ")");
    Ent("SHAPE_DEFINITION_REPRESENTATION(" + Ref(mShape) + "," + Ref(rep) + ")");
    mOut += "ENDSEC;\nEND-ISO-10303-21;\n";
    if (mPath.empty())
        return;  // built for its report only
    auto os = make_classic_ofstream(mPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + mPath);
    os << mOut;
}

}  // namespace detail
}  // namespace meshioplusplus
