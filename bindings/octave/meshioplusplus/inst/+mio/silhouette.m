function r = silhouette(mesh, varargin)
%SILHOUETTE  The outline of a part along a direction.
%   R = mio.silhouette(MESH, ...) projects the surface (or a volume mesh's
%   skin) onto the plane normal to 'normal' ([0 0 1]) and unions it exactly:
%   a chained line mesh of closed loops, holes included, ready for
%   mio.fit_curves (roadmap 7.6.3). Options: 'clearance' (0; > 0 grows the
%   outline, < 0 shrinks it), 'round_angle' (180: corners turning by less are
%   rounded), 'arc_segments' (64), 'extrude' (false: true returns the closed
%   surface of the prism over the outline), 'region' ('') and
%   'merge_tolerance' (1e-9). R is a struct: mesh, num_loops, num_outer,
%   num_holes, num_triangles, area, perimeter, min_height, max_height and
%   warnings.

    o = mio.internal.opts(varargin, {'normal', [0 0 1], 'clearance', 0, 'round_angle', 180, ...
        'arc_segments', 64, 'extrude', false, 'region', '', 'merge_tolerance', 1e-9}, ...
        'mio.silhouette');
    r = meshioplusplus_mex('mio_silhouette', mesh, double(o.normal), double(o.clearance), ...
        double(o.round_angle), double(o.arc_segments), logical(o.extrude), o.region, ...
        double(o.merge_tolerance));
    if isempty(r.warnings)
        r.warnings = cell(0, 1);
    else
        r.warnings = strsplit(r.warnings, char(10));
        r.warnings = r.warnings(:);
    end
end
