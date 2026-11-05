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
#pragma once

/**
 * @file detail/chunked_parse.hpp
 * @brief Blank-separated numbers parsed in parallel (roadmap 5.7.2). A
 * **core-private** header.
 *
 * The text is cut into chunks at blanks (`text_is_blank`, the classic-locale
 * set), so no token straddles two; the tokens of each chunk are counted in
 * parallel, the counts prefix-summed, and each chunk then parses its tokens
 * into their fixed slots in parallel. The text goes a wave of chunks at a
 * time, so a run of numbers followed by more of the file reads no further
 * than one wave past its end.
 *
 * A reader keeps its own serial parse and calls this first: `rParse(token,
 * index)` must take the whole token or refuse it (and never throw), and on
 * any refusal -- or too few tokens -- `parse_tokens` returns `npos` and the
 * reader runs its serial loop from the same position, which gives the same
 * values and the same error, at the same place, as before. Parsing a token
 * whole with the same number parser the serial loop uses yields the same
 * value bit for bit, so the parallel path changes no output.
 */

// System includes
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "parallel_workers.hpp"
#include "text_cursor.hpp"

namespace meshioplusplus {
namespace detail {

/// Below this many values a reader parses serially: the dispatch would not pay.
inline constexpr std::size_t chunked_parse_min_tokens = std::size_t{1} << 15;

/// Whether `Count` values are worth parsing in parallel on this build and now.
inline bool chunked_parse_worthwhile(std::size_t Count) {
    return Count >= chunked_parse_min_tokens && parallel_workers() > 1;
}

namespace chunked_parse_impl {

/// The bytes of one chunk: large enough that counting a chunk dwarfs the
/// dispatch, small enough that a wave overshoots a run of numbers little.
inline constexpr std::size_t kChunkBytes = std::size_t{1} << 18;

/// The first blank at or after `Pos` (or the end): a chunk boundary.
inline std::size_t next_blank(std::string_view Text, std::size_t Pos) {
    while (Pos < Text.size() && !text_is_blank(Text[Pos]))
        ++Pos;
    return Pos;
}

/// Calls `rFn(first, last)` for each token of `Text[Begin, End)`, in order,
/// until it returns false.
template <class Fn>
void for_each_token(std::string_view Text, std::size_t Begin, std::size_t End, Fn&& rFn) {
    std::size_t p = Begin;
    for (;;) {
        while (p < End && text_is_blank(Text[p]))
            ++p;
        if (p >= End)
            return;
        const std::size_t first = p;
        while (p < End && !text_is_blank(Text[p]))
            ++p;
        if (!rFn(first, p))
            return;
    }
}

}  // namespace chunked_parse_impl

/**
 * @brief Parses the `Count` tokens of `Text` that follow `Begin`, in parallel:
 * `rParse(std::string_view token, std::size_t index)` stores token `index`
 * and returns whether it took the whole token.
 * @return The position just past the last token; `npos` when a token was
 * refused or fewer than `Count` follow `Begin` (then nothing is promised about
 * what `rParse` stored, and the caller parses serially).
 */
template <class Parse>
std::size_t parse_tokens(std::string_view Text, std::size_t Begin, std::size_t Count,
                         Parse&& rParse) {
    using chunked_parse_impl::for_each_token;
    using chunked_parse_impl::kChunkBytes;
    using chunked_parse_impl::next_blank;
    if (Count == 0)
        return Begin;
    const std::size_t wave = 4 * std::max<std::size_t>(1, parallel_workers());
    std::size_t pos = Begin;
    std::size_t done = 0;  // tokens before `pos`
    std::vector<std::size_t> cuts;
    std::vector<std::size_t> first;
    while (pos < Text.size()) {
        cuts.assign(1, pos);
        while (cuts.size() <= wave && cuts.back() < Text.size())
            cuts.push_back(next_blank(Text, std::min(Text.size(), cuts.back() + kChunkBytes)));
        const std::size_t nch = cuts.size() - 1;
        first.assign(nch + 1, 0);
        parallel_for(
            nch,
            [&](std::size_t c) {
                std::size_t n = 0;
                for_each_token(Text, cuts[c], cuts[c + 1], [&](std::size_t, std::size_t) {
                    ++n;
                    return true;
                });
                first[c + 1] = n;
            },
            1);
        first[0] = done;
        for (std::size_t c = 0; c < nch; ++c)
            first[c + 1] += first[c];
        std::atomic<bool> refused{false};
        std::size_t end = std::string_view::npos;
        parallel_for(
            nch,
            [&](std::size_t c) {
                if (first[c] >= Count)
                    return;
                std::size_t index = first[c];
                for_each_token(Text, cuts[c], cuts[c + 1], [&](std::size_t b, std::size_t e) {
                    if (refused.load(std::memory_order_relaxed) ||
                        !rParse(Text.substr(b, e - b), index)) {
                        refused.store(true, std::memory_order_relaxed);
                        return false;
                    }
                    if (++index == Count) {
                        end = e;  // one chunk holds the last token
                        return false;
                    }
                    return true;
                });
            },
            1);
        if (refused.load())
            return std::string_view::npos;
        if (first[nch] >= Count)
            return end;
        done = first[nch];
        pos = cuts.back();
    }
    return std::string_view::npos;
}

/// How many tokens `Text` holds, counted in parallel chunks.
inline std::size_t count_tokens(std::string_view Text) {
    using chunked_parse_impl::for_each_token;
    using chunked_parse_impl::kChunkBytes;
    using chunked_parse_impl::next_blank;
    std::vector<std::size_t> cuts{0};
    while (cuts.back() < Text.size())
        cuts.push_back(next_blank(Text, std::min(Text.size(), cuts.back() + kChunkBytes)));
    const std::size_t nch = cuts.size() - 1;
    return parallel_reduce(
        nch, 1, std::size_t{0},
        [&](std::size_t b, std::size_t e) {
            std::size_t n = 0;
            for (std::size_t c = b; c < e; ++c)
                for_each_token(Text, cuts[c], cuts[c + 1], [&](std::size_t, std::size_t) {
                    ++n;
                    return true;
                });
            return n;
        },
        [](std::size_t a, std::size_t p) { return a + p; });
}

}  // namespace detail
}  // namespace meshioplusplus
