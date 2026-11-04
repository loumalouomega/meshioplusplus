function r = hessian(mesh, array, varargin)
%HESSIAN  Hessian (second derivative) of the scalar point data ARRAY.
%   Two mio.gradient passes: the result is (n x 9), the flattened row-major
%   3x3 Hessian, named '<array>:hessian' unless 'output' overrides it.
%   Options: 'method' ('green-gauss' or 'least-squares', forwarded to both
%   passes), 'location' ('cell' or 'point'), 'output' ('') and 'overwrite'
%   (false). Returns a struct: mesh, num_skipped and num_fallback.

    o = mio.internal.opts(varargin, {'method', 'green-gauss', 'location', 'cell', ...
        'output', '', 'overwrite', false}, 'mio.hessian');
    r = meshioplusplus_mex('mio_hessian', mesh, mio.internal.str(array, 'array'), ...
        mio.internal.str(o.method, 'method'), mio.internal.str(o.location, 'location'), ...
        mio.internal.str(o.output, 'output'), logical(o.overwrite));
end
