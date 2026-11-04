function tf = meshes_equal(a, b, varargin)
%MESHES_EQUAL  Whether two meshes are equal within a tolerance.
%   TF = mio.meshes_equal(A, B, ...) compares points, connectivity and data.
%   Options: 'atol' (0), 'rtol' (0), 'unordered' (match points by spatial
%   proximity instead of by index; false). TF is a logical scalar.

    o = mio.internal.opts(varargin, {'atol', 0, 'rtol', 0, 'unordered', false}, ...
        'mio.meshes_equal');
    tf = meshioplusplus_mex('mio_meshes_equal', a, b, double(o.atol), double(o.rtol), ...
        logical(o.unordered));
end
