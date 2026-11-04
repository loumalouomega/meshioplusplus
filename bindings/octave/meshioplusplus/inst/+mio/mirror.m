function r = mirror(mesh, varargin)
%MIRROR  Reflect a mesh across a plane.
%   R = mio.mirror(MESH, ...) reflects MESH across the plane through 'origin'
%   with 'normal', every cell right side out (roadmap 3.2.5): its nodes
%   permuted and Side region entries moved to their facets' new numbers.
%   Options:
%     'origin'              a point of the plane ([0 0 0])
%     'normal'              the plane's normal ([1 0 0])
%     'keep_original'       return the input followed by its image (false)
%     'weld'                with keep_original, share the points on the plane (true)
%     'tolerance'           on the plane: a fraction of the bounding-box diagonal (1e-9)
%     'reflect_vector_data' reflect float vector and tensor data (false)
%   R is a struct: mesh, num_flipped_cells, num_points_on_plane,
%   num_dropped_cells, num_dropped_sides and warnings (n x 1 cellstr).

    o = mio.internal.opts(varargin, {'origin', [0 0 0], 'normal', [1 0 0], ...
        'keep_original', false, 'weld', true, 'tolerance', 1e-9, ...
        'reflect_vector_data', false}, 'mio.mirror');
    r = meshioplusplus_mex('mio_mirror', mesh, double(o.origin), double(o.normal), ...
        logical(o.keep_original), logical(o.weld), double(o.tolerance), ...
        logical(o.reflect_vector_data));
    r.warnings = split_warnings(r.warnings);
end

function w = split_warnings(text)
    if isempty(text)
        w = cell(0, 1);
    else
        w = strsplit(text, char(10));
        w = w(:);
    end
end
