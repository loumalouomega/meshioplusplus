// SPDX-License-Identifier: MIT
/// @file series.cpp
/// @brief Implementation of series.hpp.

#include "series.hpp"

#include <filesystem>
#include <system_error>
#include <utility>

#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus::cli::tui {

PathSeries::PathSeries(std::vector<std::string> Paths, std::string Pattern, std::string Format) {
    mInput.mPaths = std::move(Paths);
    mInput.mPattern = std::move(Pattern);
    mInput.mFormat = std::move(Format);
    mEntries = sequence_expand(mInput);
}

Mesh PathSeries::Load(std::size_t Index) {
    return sequence_read_step(mEntries, Index, mInput.mFormat, ReadOptions{});
}

std::string PathSeries::Label(std::size_t Index) {
    if (Index >= mEntries.size())
        return "no step";
    const SequenceEntry& entry = mEntries[Index];
    std::string label = "step " + std::to_string(Index + 1) + "/" + std::to_string(mEntries.size());
    if (entry.mTimeSource != SequenceTimeSource::Index) {
        char buffer[40];
        detail::snprintf_c(buffer, sizeof(buffer), "%.6g", entry.mTime);
        label += std::string("  t=") + buffer;
    }
    return label + "  " + std::filesystem::path(entry.mPath).filename().string();
}

bool PathSeries::Refresh() {
    try {
        std::vector<SequenceEntry> fresh = sequence_expand(mInput);
        const bool changed = fresh.size() != mEntries.size();
        mEntries = std::move(fresh);
        return changed;
    } catch (const std::exception&) {
        return false;  // a glob that matches nothing for a moment: keep what we have
    }
}

std::string PathSeries::Fingerprint() {
    // A glob is expanded afresh, so a new file shows up as a new count and path.
    if (!mInput.mPattern.empty())
        Refresh();
    if (mEntries.empty())
        return "";
    const std::string& path = mEntries.back().mPath;
    std::error_code ec;
    const auto size = std::filesystem::file_size(path, ec);
    const auto time = std::filesystem::last_write_time(path, ec);
    return std::to_string(mEntries.size()) + "|" + path + "|" +
           std::to_string(ec ? 0 : static_cast<unsigned long long>(size)) + "|" +
           std::to_string(ec ? 0 : static_cast<long long>(time.time_since_epoch().count()));
}

}  // namespace meshioplusplus::cli::tui
