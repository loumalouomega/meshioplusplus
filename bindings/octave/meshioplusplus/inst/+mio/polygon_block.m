function out = polygon_block(mesh, varargin)
%POLYGON_BLOCK  A ragged polygon block (default 1) as a cell array of 1-based row
%   vectors, one per cell.

    o = mio.internal.opts(varargin, {'block', 1}, 'mio.polygon_block');
    out = meshioplusplus_mex('mio_polygon_block', mesh, mio.internal.block(o.block));
end
