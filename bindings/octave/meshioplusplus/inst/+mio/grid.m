function m = grid(dims, varargin)
%GRID  A regular hexahedral lattice built from nothing.
%   M = mio.grid(DIMS) builds a lattice of DIMS = [nx ny nz] cells. An empty
%   lattice ([0 0 0]) is a legal request. Options:
%     'origin'     lattice lo corner ([0 0 0])
%     'spacing'    cell size per axis ([1 1 1])
%     'max_cells'  refuse a grid larger than this, by name (20000000)

    o = mio.internal.opts(varargin, {'origin', [0 0 0], 'spacing', [1 1 1], ...
        'max_cells', 20000000}, 'mio.grid');
    m = meshioplusplus_mex('mio_grid', double(dims), double(o.origin), double(o.spacing), ...
        double(o.max_cells));
end
