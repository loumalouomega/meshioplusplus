function m = slice(mesh, origin, normal, varargin)
%SLICE  Cut MESH with the plane through ORIGIN with normal NORMAL.
%   Returns the cut as a new (possibly empty) mesh. Options:
%   'record_parent_ids' (false).

    o = mio.internal.opts(varargin, {'record_parent_ids', false}, 'mio.slice');
    m = meshioplusplus_mex('mio_slice', mesh, double(origin), double(normal), ...
        logical(o.record_parent_ids));
end
