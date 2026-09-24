// The refined libMesh fixtures of tests/python/meshes/libmesh (tree_quad.xd?,
// tree_hex20.xdr), written by libMesh 1.8.0 itself (LGPL-2.1; not part of
// meshio++). Build against an installed libMesh and run once:
//
//   $(libmesh-config --cxx) $(libmesh-config --cppflags --cxxflags --include) \
//       tools/libmesh_tree_driver.cpp -o lmtree $(libmesh-config --ldflags --libs)
//   ./lmtree <outdir>      # then rename quad_tree.* -> tree_quad.*, hex20_tree.xdr -> tree_hex20.xdr
//
// quad_tree: a 2x2 QUAD4 square refined uniformly, then one child again (three
// levels), a named side set, a node set, p-level 1 on every third element.
// hex_tree: a 2x1x1 HEX8 slab with the second hex refined (not committed).
// hex20_tree: a HEX20 with an edge set on a level-0 edge, refined uniformly.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include <string>
using namespace libMesh;
int main(int argc, char** argv) {
    LibMeshInit init(argc, argv);
    const std::string out = argv[1];
    {   // 2-D: a 2x2 QUAD4 square, refined once, then one child again
        Mesh mesh(init.comm());
        mesh.allow_renumbering(false);
        MeshTools::Generation::build_square(mesh, 2, 2, 0., 1., 0., 1., QUAD4);
        mesh.get_boundary_info().sideset_name(0) = "bottom";
        mesh.get_boundary_info().add_node(mesh.node_ptr(0), 7);
        mesh.get_boundary_info().nodeset_name(7) = "corner";
        mesh.subdomain_name(0) = "plate";
        MeshRefinement r(mesh);
        r.uniformly_refine(1);
        for (auto* e : mesh.active_element_ptr_range()) { e->set_refinement_flag(Elem::REFINE); break; }
        r.refine_elements();
        for (auto* e : mesh.active_element_ptr_range()) if (e->id() % 3 == 0) e->set_p_level(1);
        for (auto ext : {".xda", ".xdr"}) mesh.write(out + "/quad_tree" + ext);
    }
    {   // 3-D: a 2x1x1 HEX8 slab, the second hex refined
        Mesh mesh(init.comm());
        mesh.allow_renumbering(false);
        MeshTools::Generation::build_cube(mesh, 2, 1, 1, 0., 2., 0., 1., 0., 1., HEX8);
        mesh.get_boundary_info().sideset_name(2) = "right";
        for (auto* e : mesh.active_element_ptr_range()) if (e->id() == 1) e->set_refinement_flag(Elem::REFINE);
        MeshRefinement r(mesh);
        r.refine_elements();
        for (auto ext : {".xda", ".xdr"}) mesh.write(out + "/hex_tree" + ext);
    }
    {   // 3-D quadratic with edges: a HEX20 with an edge set, refined
        Mesh mesh(init.comm());
        mesh.allow_renumbering(false);
        MeshTools::Generation::build_cube(mesh, 1, 1, 1, 0., 1., 0., 1., 0., 1., HEX20);
        Elem* e0 = mesh.elem_ptr(0);
        mesh.get_boundary_info().add_edge(e0, 0, 11);
        mesh.get_boundary_info().sideset_name(11) = "axis";
        MeshRefinement r(mesh);
        r.uniformly_refine(1);
        for (auto ext : {".xda", ".xdr"}) mesh.write(out + "/hex20_tree" + ext);
    }
    return 0;
}
