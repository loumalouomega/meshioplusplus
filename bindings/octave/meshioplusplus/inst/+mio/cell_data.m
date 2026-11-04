function x = cell_data(mesh, name, varargin)
%CELL_DATA  A copy of cell data NAME for block BLOCK (default 1), in its stored class.

    o = mio.internal.opts(varargin, {'block', 1}, 'mio.cell_data');
    x = meshioplusplus_mex('mio_cell_data', mesh, mio.internal.str(name, 'name'), ...
        mio.internal.block(o.block));
end
