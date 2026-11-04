function r = tet_to_poly(mesh, varargin)
%TET_TO_POLY  The polyhedral (median) dual of a tetrahedral mesh.
%   R = mio.tet_to_poly(MESH, ...): one polyhedron per vertex, bounded by the
%   ring of tetrahedron and triangle centroids around each edge leaving it
%   and, on the boundary, by its share of the boundary triangles, fused while
%   their normals differ by at most 'feature_angle' degrees. A vertex on an
%   interface between tetra blocks gets one cell per block. Each tetra block
%   becomes a polyhedron block and each triangle block a polygon block; any
%   other block is dropped. Options:
%     'feature_angle'           degrees in [0, 180) (1)
%     'point_data_to_cell_data' copy each vertex's point data onto its dual
%                               cells as cell data (false)
%     'record_vertex_ids'       attach 'dual:vertex' (0-based vertex ids) (false)
%   R is a struct: mesh, vertex_map (1-based, one entry per output cell),
%   num_new_points, num_boundary_faces and num_split_vertices.

    o = mio.internal.opts(varargin, {'feature_angle', 1, 'point_data_to_cell_data', false, ...
        'record_vertex_ids', false}, 'mio.tet_to_poly');
    r = meshioplusplus_mex('mio_tet_to_poly', mesh, double(o.feature_angle), ...
        logical(o.point_data_to_cell_data), logical(o.record_vertex_ids));
end
