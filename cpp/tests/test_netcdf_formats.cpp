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
// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"

#ifdef MESHIOPLUSPLUS_HAS_NETCDF

#include "meshioplusplus/formats/exodus.hpp"

TEST(Exodus, Basic) {
    auto w = [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_exodus(p, m); };
    auto r = [](const std::string& p) { return meshioplusplus::read_exodus(p); };
    mt::roundtrip(w, r, mt::tri_mesh(), ".e");
    mt::roundtrip(w, r, mt::tet_mesh(), ".e");
    mt::roundtrip(w, r, mt::hex_mesh(), ".e");
    mt::roundtrip(w, r, mt::tri_quad_mesh(), ".e");
}

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
