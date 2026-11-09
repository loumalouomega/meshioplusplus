function r = measure(mesh, spec)
%MEASURE  Elements, dimensions and GD&T tolerances.
%   R = mio.measure(MESH, SPEC) measures MESH by the JSON spec SPEC (its text,
%   or the path of a JSON file): elements fitted from regions, seeds or
%   patches, constructed or literal; dimensions and tolerances (roadmap 7.2,
%   doc/measure.md). R is a struct: mesh (with its element table),
%   num_elements, num_dimensions, num_tolerances, num_fail, report_json and
%   warnings.

    if isempty(regexp(spec, '^\s*\{', 'once')) && exist(spec, 'file')
        spec = fileread(spec);
    end
    r = meshioplusplus_mex('mio_measure', mesh, mio.internal.str(spec, 'spec'));
    if isempty(r.warnings)
        r.warnings = cell(0, 1);
    else
        r.warnings = strsplit(r.warnings, char(10));
        r.warnings = r.warnings(:);
    end
end
