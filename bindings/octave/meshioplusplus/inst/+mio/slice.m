function m = slice(mesh, origin, normal, varargin)
%SLICE  Cut MESH with the plane through ORIGIN with normal NORMAL.
%   Returns the cut as a new (possibly empty) mesh. Options:
%   'record_parent_ids' (false), 'chain' (false: walk a surface's line
%   section into polylines with section:chain and section:closed cell data,
%   roadmap 7.6.1) and 'min_length' (0: drop shorter chains).

    o = mio.internal.opts(varargin, {'record_parent_ids', false, 'chain', false, ...
        'min_length', 0}, 'mio.slice');
    m = meshioplusplus_mex('mio_slice', mesh, double(origin), double(normal), ...
        logical(o.record_parent_ids), logical(o.chain), double(o.min_length));
end
