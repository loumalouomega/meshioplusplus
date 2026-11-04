function out = polyhedron_block(mesh, varargin)
%POLYHEDRON_BLOCK  A polyhedron block (default 1) as a cell array (one per cell) of cell
%   arrays of 1-based face rows.

    o = mio.internal.opts(varargin, {'block', 1}, 'mio.polyhedron_block');
    out = meshioplusplus_mex('mio_polyhedron_block', mesh, mio.internal.block(o.block));
end
