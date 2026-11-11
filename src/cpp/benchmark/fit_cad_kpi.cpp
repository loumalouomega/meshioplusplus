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

// fit_cad's KPI harness (doc/fit_cad_freeform_plan.md, step 1): every fixture
// of the golden, freeform and kinematic sets (and, on request, the hook of
// notebook 39) through every freeform mode and backend, one CSV row each, with
// the validation report on. tools/fit_cad_kpi.py runs it, writes
// benchmark/fit_cad_kpi.csv and the table in doc/benchmarks.md, and compares
// two runs.
//
//   meshioplusplus_fit_cad_kpi [--set golden,freeform,kinematic] [--variants]
//       [--fixtures a,b] [--modes legacy,adaptive,network,auto]
//       [--backends native,occt] [--hook PATH] [--option key=value]...
//       [--tag TEXT] [--step-dir DIR] [--csv FILE]
//
// An `--option` sets a FitCadOptions member for every run (the plan's new
// options, so two runs differ only in it); `--tag` names the run in the CSV.

// System includes
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/skin.hpp"

// Test fixtures (tests/cpp, on this target's include path)
#include "fit_cad_fixtures.hpp"

namespace {

using meshioplusplus::FitCadBackend;
using meshioplusplus::FitCadFreeformMode;
using meshioplusplus::FitCadOptions;
using meshioplusplus::FitCadResult;
using meshioplusplus::Mesh;

struct KpiFixture {
    std::string mName;
    std::string mSet;
    Mesh mMesh;
    double mTolerance = 0.0;  // 0: fit_cad's own
};

std::vector<std::string> kpi_split(const std::string& rText) {
    std::vector<std::string> out;
    auto ss = meshioplusplus::detail::make_classic_stringstream(rText);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty())
            out.push_back(item);
    return out;
}

bool kpi_has(const std::vector<std::string>& rList, const std::string& rItem) {
    for (const std::string& s : rList)
        if (s == rItem)
            return true;
    return false;
}

// The hook of notebook 39: the tetrahedral deck's skin, segmented by the
// primitive method and fitted at 1e-3 of its diagonal, as the notebook does.
KpiFixture kpi_hook(const std::string& rPath) {
    meshioplusplus::ReadOptions ro;
    const Mesh deck = meshioplusplus::registry_read(rPath, "", ro);
    const Mesh skin = meshioplusplus::extract_skin(deck);
    const double tol = 1e-3 * fcfix::diagonal(skin);
    meshioplusplus::SegmentOptions so;
    so.mMethod = meshioplusplus::SegmentMethod::Primitive;
    meshioplusplus::FitPrimitivesOptions fo;
    fo.mTolerance = tol;
    fo.mRelative = false;
    Mesh fitted = meshioplusplus::fit_primitives(meshioplusplus::segment(skin, so).mMesh, fo).mMesh;
    return {"hook", "hook", std::move(fitted), 0.0};
}

// Sets a FitCadOptions member from `key=value`; false for an unknown key.
bool kpi_option(FitCadOptions& rO, const std::string& rKey, const std::string& rValue) {
    const bool on = rValue == "1" || rValue == "true" || rValue == "on";
    static const std::map<std::string,
                          std::function<void(FitCadOptions&, const std::string&, bool)>>
        kSetters = {
            {"validation_acceptance",
             [](FitCadOptions& o, const std::string&, bool b) { o.mValidationAcceptance = b; }},
            {"constrain_primitives",
             [](FitCadOptions& o, const std::string&, bool b) { o.mConstrainPrimitives = b; }},
            {"refine_segmentation",
             [](FitCadOptions& o, const std::string&, bool b) { o.mRefineSegmentation = b; }},
            {"local_refine",
             [](FitCadOptions& o, const std::string& v, bool) {
                 o.mLocalRefine = v == "thb" ? meshioplusplus::FitCadLocalRefine::Thb
                                             : meshioplusplus::FitCadLocalRefine::Knots;
             }},
            {"freeform_fit_schedule",
             [](FitCadOptions& o, const std::string& v, bool) {
                 o.mFreeformFitSchedule =
                     v == "continuation" ? meshioplusplus::FitCadFreeformFitSchedule::Continuation
                                         : meshioplusplus::FitCadFreeformFitSchedule::Legacy;
             }},
            {"freeform_objective",
             [](FitCadOptions& o, const std::string& v, bool) {
                 o.mFreeformObjective = v == "pdm"   ? meshioplusplus::FitCadFreeformObjective::Pdm
                                        : v == "tdm" ? meshioplusplus::FitCadFreeformObjective::Tdm
                                        : v == "sdm"
                                            ? meshioplusplus::FitCadFreeformObjective::Sdm
                                            : meshioplusplus::FitCadFreeformObjective::Default;
             }},
            {"resegment_rounds", [](FitCadOptions& o, const std::string& v,
                                    bool) { o.mResegmentRounds = std::stoll(v); }},
        };
    const auto it = kSetters.find(rKey);
    if (it == kSetters.end())
        return false;
    it->second(rO, rValue, on);
    return true;
}

const char* kpi_mode_name(FitCadFreeformMode M) {
    static const char* const kNames[] = {"legacy", "adaptive", "network", "auto"};
    return kNames[static_cast<int>(M) & 3];
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> sets{"golden", "freeform", "kinematic"}, only, modes, backends;
    std::vector<std::pair<std::string, std::string>> options;
    std::string hook, csv, tag = "default", step_dir;
    bool variants = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const auto next = [&]() -> std::string {
            if (i + 1 >= argc)
                throw std::invalid_argument(a + " needs a value");
            return argv[++i];
        };
        if (a == "--set")
            sets = kpi_split(next());
        else if (a == "--fixtures")
            only = kpi_split(next());
        else if (a == "--modes")
            modes = kpi_split(next());
        else if (a == "--backends")
            backends = kpi_split(next());
        else if (a == "--variants")
            variants = true;
        else if (a == "--hook")
            hook = next();
        else if (a == "--csv")
            csv = next();
        else if (a == "--tag")
            tag = next();
        else if (a == "--step-dir")
            step_dir = next();
        else if (a == "--option") {
            const std::string kv = next();
            const std::size_t eq = kv.find('=');
            if (eq == std::string::npos)
                throw std::invalid_argument("--option takes key=value");
            options.emplace_back(kv.substr(0, eq), kv.substr(eq + 1));
        } else if (a == "--help" || a == "-h") {
            std::cout << "see the file comment of src/cpp/benchmark/fit_cad_kpi.cpp\n";
            return 0;
        } else {
            throw std::invalid_argument("unknown argument " + a);
        }
    }
    if (modes.empty())
        modes = {"legacy", "adaptive", "network", "auto"};
    if (backends.empty())
        backends = meshioplusplus::fit_cad_has_occt() ? std::vector<std::string>{"native", "occt"}
                                                      : std::vector<std::string>{"native"};
    if (step_dir.empty())
        step_dir = (std::filesystem::temp_directory_path() / "meshioplusplus_fit_cad_kpi").string();
    std::filesystem::create_directories(step_dir);

    // The fixtures.
    std::vector<KpiFixture> fixtures;
    const auto want = [&](const std::string& rName) {
        return only.empty() || kpi_has(only, rName);
    };
    if (kpi_has(sets, "golden"))
        for (fcfix::Fixture& f : fcfix::golden_fixtures())
            if (want(f.mName))
                fixtures.push_back({f.mName, "golden", std::move(f.mMesh), 0.0});
    if (kpi_has(sets, "freeform"))
        for (fcfix::FreeformFixture& f : fcfix::freeform_fixtures(variants))
            if (want(f.mName))
                fixtures.push_back({f.mName, "freeform", std::move(f.mMesh), f.mTolerance});
    if (kpi_has(sets, "kinematic"))
        for (fcfix::FreeformFixture& f : fcfix::kinematic_fixtures())
            if (want(f.mName))
                fixtures.push_back({f.mName, "kinematic", std::move(f.mMesh), f.mTolerance});
    if (!hook.empty())
        fixtures.push_back(kpi_hook(hook));

    auto file = meshioplusplus::detail::make_classic_ofstream();
    if (!csv.empty())
        file.open(csv);
    std::ostream& out = csv.empty() ? std::cout : file;
    out << "tag,set,fixture,backend,mode,valid,closed,brep_checked,tolerance,volume_ratio,"
           "facet_area_fraction,num_faces,num_plane,num_sphere,num_cylinder,num_cone,num_torus,"
           "num_bspline,num_fillet,num_helical,num_faceted_faces,num_demoted,"
           "num_network_patches,total_ctrl,mesh_to_cad_max,mesh_to_cad_p95,cad_to_mesh_max,"
           "cad_to_mesh_p95,normal_max_deg,normal_p95_deg,max_g1_deg,step_round_trip,"
           "step_volume_ratio,seconds\n";
    for (const KpiFixture& f : fixtures)
        for (const std::string& backend : backends)
            for (const std::string& mode : modes) {
                FitCadOptions o;
                o.mTolerance = f.mTolerance;
                o.mBackend = backend == "occt" ? FitCadBackend::Occt : FitCadBackend::Native;
                o.mFreeformMode = mode == "adaptive"  ? FitCadFreeformMode::Adaptive
                                  : mode == "network" ? FitCadFreeformMode::Network
                                  : mode == "auto"    ? FitCadFreeformMode::Auto
                                                      : FitCadFreeformMode::Legacy;
                o.mValidationReport = true;
                for (const auto& [k, v] : options)
                    if (!kpi_option(o, k, v))
                        throw std::invalid_argument("unknown --option " + k);
                const std::string step = (std::filesystem::path(step_dir) /
                                          (f.mName + "_" + backend + "_" + mode + ".step"))
                                             .string();
                const auto t0 = std::chrono::steady_clock::now();
                FitCadResult r;
                try {
                    r = meshioplusplus::fit_cad(f.mMesh, step, o);
                } catch (const std::exception& rErr) {
                    std::cerr << f.mName << " " << backend << " " << mode
                              << ": failed: " << rErr.what() << "\n";
                    continue;
                }
                const double secs =
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                char line[2048];
                std::snprintf(
                    line, sizeof line,
                    "%s,%s,%s,%s,%s,%d,%d,%d,%.6g,%.6g,%.6g,%lld,%lld,%lld,%lld,%lld,%lld,%lld,"
                    "%lld,%lld,%lld,%lld,%lld,%lld,%.6g,%.6g,%.6g,%.6g,%.4g,%.4g,%.4g,%d,%.6g,"
                    "%.3f\n",
                    tag.c_str(), f.mSet.c_str(), f.mName.c_str(), backend.c_str(),
                    kpi_mode_name(o.mFreeformMode), r.mValid ? 1 : 0, r.mClosed ? 1 : 0,
                    r.mBrepChecked ? 1 : 0, r.mTolerance, r.mVolumeRatio, r.mFacetAreaFraction,
                    static_cast<long long>(r.mCounts.mNumFaces),
                    static_cast<long long>(r.mCounts.mNumPlane),
                    static_cast<long long>(r.mCounts.mNumSphere),
                    static_cast<long long>(r.mCounts.mNumCylinder),
                    static_cast<long long>(r.mCounts.mNumCone),
                    static_cast<long long>(r.mCounts.mNumTorus),
                    static_cast<long long>(r.mNumBspline), static_cast<long long>(r.mNumFillet),
                    static_cast<long long>(r.mNumHelical),
                    static_cast<long long>(r.mCounts.mNumFacetedFaces),
                    static_cast<long long>(r.mCounts.mNumDemoted),
                    static_cast<long long>(r.mNumNetworkPatches),
                    static_cast<long long>(r.mTotalFreeformControlPoints), r.mMeshToCadMax,
                    r.mMeshToCadP95, r.mCadToMeshMax, r.mCadToMeshP95, r.mNormalDeviationMaxDeg,
                    r.mNormalDeviationP95Deg, r.mMaxG1AngleDeg,
                    r.mStepRoundTripChecked ? (r.mStepRoundTripValid ? 1 : 0) : -1,
                    r.mStepRoundTripVolumeRatio, secs);
                out << line;
                out.flush();
                std::cerr << f.mName << " " << backend << " " << mode << ": facets "
                          << r.mFacetAreaFraction << ", " << secs << " s\n";
            }
    return 0;
}
