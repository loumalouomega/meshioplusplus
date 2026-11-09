function r = flip_normals(mesh, varargin)
%FLIP_NORMALS  Reverse the winding of cells.
%   R = mio.flip_normals(MESH, ...) reverses every cell, or the cells of the
%   cell region 'region' (''), with no point moved; Side entries move with
%   their cells (roadmap 7.8.2). R is a struct: mesh and num_flipped.

    o = mio.internal.opts(varargin, {'region', ''}, 'mio.flip_normals');
    r = meshioplusplus_mex('mio_flip_normals', mesh, o.region);
end
