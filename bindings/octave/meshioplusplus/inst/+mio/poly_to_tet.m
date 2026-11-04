function r = poly_to_tet(mesh, varargin)
%POLY_TO_TET  Decompose a polyhedral mesh into tetrahedra.
%   R = mio.poly_to_tet(MESH, ...) works like mio.convert_cells(MESH,
%   'simplexify'), except that a polyhedron that already is a tetrahedron
%   becomes that one tetra. Option 'record_parent_ids' (attach
%   'convert:parent_cell'; false). R is a struct: mesh, point_map, cell_maps
%   (as mio.convert_cells), num_exact_tets and num_fanned_cells.

    o = mio.internal.opts(varargin, {'record_parent_ids', false}, 'mio.poly_to_tet');
    r = meshioplusplus_mex('mio_poly_to_tet', mesh, logical(o.record_parent_ids));
end
