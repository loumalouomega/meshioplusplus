function m = crop_plane(mesh, point, normal, varargin)
%CROP_PLANE  Keep the cells on the NORMAL side of a plane; returns a new mesh.
%   The plane passes through POINT (3-vector) with normal NORMAL (3-vector).
%   Options: 'mode' ('all' or 'any') and 'record_ids' (false).

    o = mio.internal.opts(varargin, {'mode', 'all', 'record_ids', false}, 'mio.crop_plane');
    m = meshioplusplus_mex('mio_crop_plane', mesh, double(point), double(normal), ...
        mio.internal.choice(o.mode, {'all', 'any'}, 'mode'), logical(o.record_ids));
end
