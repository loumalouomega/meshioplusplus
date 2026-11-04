function r = tessellate(mesh, varargin)
%TESSELLATE  Tessellate the curved cells of a mesh.
%   R = mio.tessellate(MESH, ...) samples every triangle6, quad8, quad9,
%   tetra10 or hexahedron27 cell on a reference lattice mapped through its
%   own quadratic shape functions; curved surface cells become one
%   'triangle' block and curved volume cells one 'tetra' block, after the
%   pass-through blocks. The output is conforming. See doc/tessellation.md.
%   Options:
%     'levels'         divisions per axis for quad8/quad9/hexahedron27, or
%                      recursive splits for triangle6/tetra10 (2^levels per
%                      edge); >= 0 (2)
%     'curved'         false passes every cell through unchanged (true)
%     'fields'         interpolate point data and copy cell data (true)
%     'record_stencil' also store the (num_points x 27) point data
%                      'tessellate:stencil' and 'tessellate:weights' (false)
%   R is a struct: mesh (with point data 'tessellate:source_point', the
%   0-based source point or -1, and cell data 'tessellate:source_cell' and
%   'tessellate:sub_index'), num_source_points, num_source_cells,
%   num_curved_source_cells, num_pass_through_source_cells, num_points,
%   num_cells, num_boundary_facets, num_interior_facets and
%   num_facets_with_bad_count (zero for a conforming tessellation).

    o = mio.internal.opts(varargin, {'levels', 2, 'curved', true, 'fields', true, ...
        'record_stencil', false}, 'mio.tessellate');
    r = meshioplusplus_mex('mio_tessellate', mesh, double(o.levels), logical(o.curved), ...
        logical(o.fields), logical(o.record_stencil));
end
