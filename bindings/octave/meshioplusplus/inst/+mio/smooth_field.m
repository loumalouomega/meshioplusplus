function r = smooth_field(mesh, array, varargin)
%SMOOTH_FIELD  Screened-Poisson smoothing of a point-data field.
%   Options: 'output' (''), 'length' (0.0), 'tolerance' (1e-12),
%   'max_iterations' (0).

    o = mio.internal.opts(varargin, {'output', '', 'length', 0.0, ...
        'tolerance', 1e-12, 'max_iterations', 0}, 'mio.smooth_field');
    r = meshioplusplus_mex('mio_smooth_field', mesh, char(array), ...
        char(o.output), double(o.length), double(o.tolerance), ...
        double(o.max_iterations));
end
