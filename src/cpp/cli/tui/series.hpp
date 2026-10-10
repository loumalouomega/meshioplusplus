// SPDX-License-Identifier: MIT
/// @file series.hpp
/// @brief A time series of meshes for the interactive viewer, one alive at a time.
#ifndef MESHIOPLUSPLUS_CLI_TUI_SERIES_HPP
#define MESHIOPLUSPLUS_CLI_TUI_SERIES_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/sequence.hpp"

namespace meshioplusplus::cli::tui {

/// What the viewer steps through: a count, one mesh at a time, a label for the
/// status line, and a way to notice that a live run wrote another step.
class TuiSeries {
public:
    virtual ~TuiSeries() = default;
    virtual std::size_t Count() = 0;
    /// Read one step. Throws (`ReadError` for a half-written file) rather than
    /// returning a partial mesh.
    virtual Mesh Load(std::size_t Index) = 0;
    /// `step 3/12  t=0.25  out_0003.vtu`.
    virtual std::string Label(std::size_t Index) = 0;
    /// Look again for steps; true when the count changed.
    virtual bool Refresh() = 0;
    /// Changes whenever the newest step's file does (path, size, modification
    /// time) or a step appears: the viewer waits for it to stop changing before
    /// it reads.
    virtual std::string Fingerprint() = 0;
};

/// A series from explicit paths or a glob (`out_*.vtu`, natural-numeric order),
/// through the library's sequence machinery.
class PathSeries : public TuiSeries {
public:
    /// Exactly one of `rPaths` and `rPattern` is non-empty. Throws what
    /// `sequence_expand` throws (nothing matched, a missing file).
    PathSeries(std::vector<std::string> Paths, std::string Pattern, std::string Format);

    std::size_t Count() override { return mEntries.size(); }
    Mesh Load(std::size_t Index) override;
    std::string Label(std::size_t Index) override;
    bool Refresh() override;
    std::string Fingerprint() override;

private:
    SequenceInput mInput;
    std::vector<SequenceEntry> mEntries;
};

}  // namespace meshioplusplus::cli::tui

#endif  // MESHIOPLUSPLUS_CLI_TUI_SERIES_HPP
