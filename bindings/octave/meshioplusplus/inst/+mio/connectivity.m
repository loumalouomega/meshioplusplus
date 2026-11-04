function out = connectivity(mesh, varargin)
%CONNECTIVITY  Connectivity of block BLOCK (default 1): num_cells x nodes_per_cell,
%   1-based node indices, one cell per row.

    o = mio.internal.opts(varargin, {'block', 1}, 'mio.connectivity');
    out = meshioplusplus_mex('mio_connectivity', mesh, mio.internal.block(o.block));
end
