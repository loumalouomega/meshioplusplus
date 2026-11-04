function r = tetrahedralize(mesh, varargin)
%TETRAHEDRALIZE  Fill a closed triangle surface with tetrahedra.
%   R = mio.tetrahedralize(MESH, ...) fills the closed, manifold,
%   consistently oriented triangle surface MESH with tetrahedra (roadmap
%   6.7): through Netgen its triangles and points are kept exactly, through
%   the native octree backend its boundary is fitted within
%   hausdorff_to_surface. Each shell bounds a region, meshed on its own. The
%   result has a 'tetra' block, a 'triangle' block and the cell data
%   'tetrahedralize:region'. 'auto' takes Netgen when mio.has_netgen is true
%   and the octree backend otherwise; 'tetgen' and 'gmsh' run those
%   installed programs (see mio.tetrahedralize_backends). Options:
%     'backend'        'auto', 'netgen', 'octree', or another mesher's name ('auto')
%     'max_size'       the largest edge asked of the tetrahedra; 0 is none (0)
%     'max_volume'     the largest volume; 0 is none (0)
%     'grading'        how fast the size grows away from the surface, in [0, 1] (0.3)
%     'optimize_steps' passes of the volume optimisation, 0 to 100 (3)
%     'max_cells'      refuse a result of more tetrahedra (5e7)
%     'face_array'     the cell data naming each triangle's CAD face ('cad:face')
%   R is a struct: mesh, num_tetrahedra, num_points, num_new_points,
%   num_regions, num_boundary_triangles, boundary_preserved,
%   hausdorff_to_surface, min_dihedral_deg, max_dihedral_deg, volume and
%   backend.

    o = mio.internal.opts(varargin, {'backend', 'auto', 'max_size', 0, 'max_volume', 0, ...
        'grading', 0.3, 'optimize_steps', 3, 'max_cells', 5e7, 'face_array', 'cad:face'}, ...
        'mio.tetrahedralize');
    r = meshioplusplus_mex('mio_tetrahedralize', mesh, mio.internal.str(o.backend, 'backend'), ...
        double(o.max_size), double(o.max_volume), double(o.grading), ...
        double(o.optimize_steps), double(o.max_cells), ...
        mio.internal.str(o.face_array, 'face_array'));
end
