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

// check_quality. See operations/quality_gate.hpp for the contract.

// System includes
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/quality_gate.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kQgPrefix = "meshio++: quality gate: ";

// compute_quality's metrics, in its fixed order (pinned against the report
// by tests/cpp/test_quality_gate.cpp).
const char* const kQgMetrics[] = {
    "quality:volume",   "quality:scaled_jacobian", "quality:aspect_ratio",
    "quality:skewness", "quality:min_angle",       "quality:max_angle",
    "quality:warpage",  "quality:min_dihedral",    "quality:max_dihedral",
    "quality:inverted", "quality:degenerate",
};

std::string qg_full_name(const std::string& rMetric) {
    const std::string full = rMetric.rfind("quality:", 0) == 0 ? rMetric : "quality:" + rMetric;
    for (const char* m : kQgMetrics)
        if (full == m)
            return full;
    std::string names;
    for (const char* m : kQgMetrics)
        names += (names.empty() ? "" : ", ") + std::string(m + 8);
    throw std::invalid_argument(std::string(kQgPrefix) + "unknown metric '" + rMetric +
                                "' (expected one of " + names + ")");
}

std::string qg_trim(const std::string& rS) {
    std::size_t b = 0;
    std::size_t e = rS.size();
    while (b < e && std::isspace(static_cast<unsigned char>(rS[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(rS[e - 1])))
        --e;
    return rS.substr(b, e - b);
}

double qg_number(const std::string& rText, const std::string& rClause) {
    const std::string t = qg_trim(rText);
    const char* end = nullptr;
    const double v = detail::parse_double(t.c_str(), end);
    if (t.empty() || end != t.c_str() + t.size() || !std::isfinite(v))
        throw std::invalid_argument(std::string(kQgPrefix) + "'" + rClause +
                                    "': expected a number, got '" + t + "'");
    return v;
}

std::string qg_g(double v) {
    char buf[32];
    detail::snprintf_c(buf, sizeof(buf), "%g", v);
    return buf;
}

std::string qg_name(const std::string& rMetric, double Min, double Max, double Fraction) {
    const std::string m = rMetric.substr(8);
    std::string s;
    if (!std::isnan(Min) && !std::isnan(Max))
        s = qg_g(Min) + " <= " + m + " <= " + qg_g(Max);
    else if (!std::isnan(Min))
        s = m + " >= " + qg_g(Min);
    else
        s = m + " <= " + qg_g(Max);
    if (Fraction > 0.0)
        s += " @ " + qg_g(Fraction * 100.0) + "%";
    return s;
}

}  // namespace

std::vector<QualityThreshold> parse_quality_thresholds(const std::string& rText) {
    // Strip `#` comments to the end of their line; a newline separates
    // clauses like `;` and `,` do, so a gate file is one clause per line.
    std::string rSpec;
    bool comment = false;
    for (const char ch : rText) {
        if (ch == '\n') {
            comment = false;
            rSpec += ';';
        } else if (ch == '#') {
            comment = true;
        } else if (!comment) {
            rSpec += ch;
        }
    }
    std::vector<QualityThreshold> out;
    std::size_t start = 0;
    while (start <= rSpec.size()) {
        std::size_t stop = rSpec.find_first_of(";,", start);
        if (stop == std::string::npos)
            stop = rSpec.size();
        const std::string clause = qg_trim(rSpec.substr(start, stop - start));
        start = stop + 1;
        if (clause.empty())
            continue;
        const std::size_t ge = clause.find(">=");
        const std::size_t le = clause.find("<=");
        if ((ge == std::string::npos) == (le == std::string::npos))
            throw std::invalid_argument(std::string(kQgPrefix) + "'" + clause +
                                        "': expected METRIC >= VALUE or METRIC <= VALUE");
        const std::size_t op = ge != std::string::npos ? ge : le;
        QualityThreshold t;
        t.mMetric = qg_full_name(qg_trim(clause.substr(0, op)));
        std::string rest = clause.substr(op + 2);
        const std::size_t at = rest.find('@');
        if (at != std::string::npos) {
            std::string frac = qg_trim(rest.substr(at + 1));
            rest = rest.substr(0, at);
            const bool percent = !frac.empty() && frac.back() == '%';
            if (percent)
                frac.pop_back();
            t.mMaxFraction = qg_number(frac, clause) / (percent ? 100.0 : 1.0);
        }
        const double v = qg_number(rest, clause);
        (ge != std::string::npos ? t.mMin : t.mMax) = v;
        out.push_back(t);
    }
    return out;
}

QualityGateResult check_quality(const Mesh& rMesh, const QualityGateOptions& rOptions) {
    QualityGateResult result;
    result.mReport = compute_quality(rMesh);
    const QualityReport& rep = result.mReport;

    for (const QualityThreshold& t : rOptions.mThresholds) {
        QualityCheck c;
        c.mMetric = qg_full_name(t.mMetric);
        c.mMin = t.mMin;
        c.mMax = t.mMax;
        c.mMaxFraction = t.mMaxFraction;
        if (std::isnan(t.mMin) && std::isnan(t.mMax))
            throw std::invalid_argument(std::string(kQgPrefix) + "threshold on '" + t.mMetric +
                                        "' has neither a minimum nor a maximum");
        if (!(t.mMaxFraction >= 0.0 && t.mMaxFraction <= 1.0))
            throw std::invalid_argument(std::string(kQgPrefix) + "threshold on '" + t.mMetric +
                                        "': the allowed fraction must lie in [0, 1]");
        c.mName = qg_name(c.mMetric, t.mMin, t.mMax, t.mMaxFraction);

        const std::vector<NDArray>* arrays = nullptr;
        for (const auto& kv : rep.mCellArrays)
            if (kv.first == c.mMetric)
                arrays = &kv.second;
        // The margin by which a value clears the bound (negative: violates it).
        auto margin = [&](double v) {
            double m = std::numeric_limits<double>::infinity();
            if (!std::isnan(t.mMin))
                m = std::min(m, v - t.mMin);
            if (!std::isnan(t.mMax))
                m = std::min(m, t.mMax - v);
            return m;
        };
        double worst_margin = std::numeric_limits<double>::infinity();
        std::int64_t cell = 0;
        if (arrays != nullptr) {
            for (const NDArray& a : *arrays) {
                const double* d = a.As<double>();
                for (std::size_t i = 0; i < a.Size(); ++i, ++cell) {
                    const double v = d[i];
                    if (!std::isfinite(v))
                        continue;
                    ++c.mEvaluated;
                    const double m = margin(v);
                    if (m < 0.0)
                        ++c.mViolations;
                    if (m < worst_margin) {
                        worst_margin = m;
                        c.mWorst = v;
                        c.mWorstCell = cell;
                    }
                }
            }
        }
        c.mFraction = c.mEvaluated > 0
                          ? static_cast<double>(c.mViolations) / static_cast<double>(c.mEvaluated)
                          : 0.0;
        c.mPassed = c.mViolations == 0 || c.mFraction <= t.mMaxFraction;
        if (c.mEvaluated == 0)
            log::warn(
                "quality gate: '{}' applies to no cell of this mesh (the metric is not defined "
                "for its cell types), so it passes vacuously",
                c.mName);
        result.mPassed = result.mPassed && c.mPassed;
        result.mChecks.push_back(c);
    }

    auto count_check = [&](const char* pWhat, std::int64_t Count, std::int64_t Limit) {
        if (Limit < 0)
            return;
        QualityCheck c;
        c.mName = std::string(pWhat) + " <= " + std::to_string(Limit);
        c.mMetric = pWhat;
        c.mMax = static_cast<double>(Limit);
        c.mEvaluated = rep.mNumCells;
        c.mViolations = Count;
        c.mFraction = rep.mNumCells > 0
                          ? static_cast<double>(Count) / static_cast<double>(rep.mNumCells)
                          : 0.0;
        c.mWorst = static_cast<double>(Count);
        c.mPassed = Count <= Limit;
        result.mPassed = result.mPassed && c.mPassed;
        result.mChecks.push_back(c);
    };
    count_check("inverted", rep.mNumInverted, rOptions.mMaxInverted);
    count_check("degenerate", rep.mNumDegenerate, rOptions.mMaxDegenerate);
    return result;
}

std::string quality_gate_summary(const QualityGateResult& rResult) {
    std::string s = std::string("quality gate: ") + (rResult.mPassed ? "PASS" : "FAIL") + " (" +
                    std::to_string(rResult.mReport.mNumCells) + " cells)\n";
    for (const QualityCheck& c : rResult.mChecks) {
        s += std::string("  ") + (c.mPassed ? "pass" : "FAIL") + "  " + c.mName + ": " +
             std::to_string(c.mViolations) + " of " + std::to_string(c.mEvaluated) +
             " cells violate";
        if (c.mMetric.rfind("quality:", 0) == 0 && c.mWorstCell >= 0)
            s += ", worst " + qg_g(c.mWorst) + " at cell " + std::to_string(c.mWorstCell);
        s += "\n";
    }
    return s;
}

}  // namespace meshioplusplus
