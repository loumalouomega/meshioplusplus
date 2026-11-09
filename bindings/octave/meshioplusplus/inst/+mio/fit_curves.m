function r = fit_curves(sections, varargin)
%FIT_CURVES  Lines, arcs and circles from chained sections.
%   R = mio.fit_curves(SECTIONS, ...) splits each chain of a chained line mesh
%   (mio.slice(..., 'chain', true) or mio.sections) into the runs a line or a
%   circular arc explains within the tolerance, in its plane's 2-D frame
%   (roadmap 7.6.2). Options: 'tolerance' (0: 1e-3 of the diagonal),
%   'snap_angle' (0), 'normal' ([0 0 0]: from the points), 'min_arc_points'
%   (4) and 'min_arc_sweep' (5 degrees). R is a struct: mesh, pieces (an n x 1
%   struct array: type, start, end, centre, radius, sweep, residual,
%   num_points, chain, frame), num_lines, num_arcs, num_circles,
%   max_residual, tolerance and warnings.

    o = mio.internal.opts(varargin, {'tolerance', 0, 'snap_angle', 0, 'normal', [0 0 0], ...
        'min_arc_points', 4, 'min_arc_sweep', 5}, 'mio.fit_curves');
    r = meshioplusplus_mex('mio_fit_curves', sections, double(o.tolerance), ...
        double(o.snap_angle), double(o.normal), double(o.min_arc_points), ...
        double(o.min_arc_sweep));
    if isempty(r.warnings)
        r.warnings = cell(0, 1);
    else
        r.warnings = strsplit(r.warnings, char(10));
        r.warnings = r.warnings(:);
    end
end
