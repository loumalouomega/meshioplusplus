function r = subdivide(mesh, varargin)
%SUBDIVIDE  Polyhedral subdivision: each cell split into one child per face.
%   R = mio.subdivide(MESH, ...) with option 'record_parent_ids' (attach
%   'subdivide:parent_cell'; false). R is a struct: mesh and cell_maps (one
%   1-based n x 1 first-child map per input block). There is no point map:
%   subdivide never prunes or renumbers a point.

    o = mio.internal.opts(varargin, {'record_parent_ids', false}, 'mio.subdivide');
    r = meshioplusplus_mex('mio_subdivide', mesh, logical(o.record_parent_ids));
end
