// License: MIT License, meshio++ default license: LICENSE
#pragma once

// Core-private view counterparts. The installed owning API stays unchanged.
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "meshioplusplus/detail/keyword_card.hpp"

namespace meshioplusplus::detail {

std::vector<std::string_view> split_card_view(std::string_view Line,
                                              const std::vector<CardField>& rLayout, CardMode Mode);
std::vector<std::string_view> split_fixed_view(std::string_view Line,
                                               const std::vector<CardField>& rFields);
std::int64_t card_to_int_view(std::string_view Text, const std::string& rWhere,
                              const std::string& rFormat = "LS-DYNA");
double card_to_real_view(std::string_view Text, const std::string& rWhere,
                         const std::string& rFormat = "LS-DYNA");

}  // namespace meshioplusplus::detail
