function r = reorder(mesh, varargin)
%REORDER  Renumber points and cells for locality.
%   R = mio.reorder(MESH, METHOD) with METHOD 'rcm' (default), 'morton' or
%   'hilbert'. R is a struct: mesh (the renumbered mesh), node_perm (n x 1,
%   1-based) and cell_perms (a cell array with one 1-based n x 1 permutation
%   per cell block).

    o = mio.internal.opts(varargin, {'method', 'rcm'}, 'mio.reorder');
    r = meshioplusplus_mex('mio_reorder', mesh, mio.internal.str(o.method, 'method'));
end
