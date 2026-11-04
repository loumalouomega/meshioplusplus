function m = transform(mesh, matrix, varargin)
%TRANSFORM  Apply a 4x4 affine MATRIX to MESH's points; returns a new mesh.
%   A point p maps to MATRIX * [p; 1]. Topology is untouched.
%   Options: 'rotate_vector_data' (false) also rotates vector (trailing dim 3)
%   and tensor (trailing dim 9) point data by the transform's 3x3 linear
%   block.

    o = mio.internal.opts(varargin, {'rotate_vector_data', false}, 'mio.transform');
    m = meshioplusplus_mex('mio_transform', mesh, double(matrix), logical(o.rotate_vector_data));
end
