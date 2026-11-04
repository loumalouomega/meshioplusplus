function r = reconstruct_surface(mesh, varargin)
%RECONSTRUCT_SURFACE  A triangle surface through an oriented point cloud.
%   Options: 'normals' ('normals' point data; estimated when absent unless
%   'estimate_normals' is false), 'k' (8), 'cell_size' (0), 'resolution'
%   (0), 'max_distance' (0), 'estimate_normals' (true), 'normals_k' (16) and
%   'max_cells' (2e7). Returns a struct: mesh, num_samples, spacing,
%   cell_size, max_distance, num_lattice_points, num_undefined,
%   num_triangles, num_components and estimated_normals.

    o = mio.internal.opts(varargin, {'normals', 'normals', 'k', 8, 'cell_size', 0, ...
        'resolution', 0, 'max_distance', 0, 'estimate_normals', true, 'normals_k', 16, ...
        'max_cells', 2e7}, 'mio.reconstruct_surface');
    r = meshioplusplus_mex('mio_reconstruct_surface', mesh, ...
        mio.internal.str(o.normals, 'normals'), double(o.k), double(o.cell_size), ...
        double(o.resolution), double(o.max_distance), logical(o.estimate_normals), ...
        double(o.normals_k), double(o.max_cells));
end
