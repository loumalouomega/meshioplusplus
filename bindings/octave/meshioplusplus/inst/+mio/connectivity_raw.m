function out = connectivity_raw(mesh, varargin)
%CONNECTIVITY_RAW  Like mio.connectivity, but 0-based (the C layout).

    o = mio.internal.opts(varargin, {'block', 1}, 'mio.connectivity_raw');
    out = meshioplusplus_mex('mio_connectivity_raw', mesh, mio.internal.block(o.block));
end
