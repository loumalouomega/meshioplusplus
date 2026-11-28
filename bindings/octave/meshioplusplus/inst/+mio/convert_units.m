function r = convert_units(mesh, varargin)
%CONVERT_UNITS  Convert declared fields to a consistent unit system.
%   R = mio.convert_units(MESH, ...) converts every declared field and the
%   points to 'system' ('SI', 'mm-t-s', 'mm-kg-ms', 'cm-g-s', 'mm-g-ms');
%   'strict' (true) fails on an undeclared array. Roadmap 9.10.1.
%   Returns a struct: mesh, num_arrays_converted, warnings.

    o = mio.internal.opts(varargin, {'system', 'SI', 'strict', true}, 'mio.convert_units');
    r = meshioplusplus_mex('mio_convert_units', mesh, mio.internal.str(o.system, 'system'), ...
        logical(o.strict));
    if isempty(r.warnings)
        r.warnings = {};
    else
        r.warnings = strsplit(r.warnings, char(10));
    end
end
