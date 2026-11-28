function r = subdivide_surface(mesh, varargin)
%SUBDIVIDE_SURFACE  Loop or Catmull-Clark subdivision of a surface.
%   R = mio.subdivide_surface(MESH, ...) subdivides with 'scheme' ('loop' or
%   'catmull_clark'), 'levels' (1), 'crease_angle' (60 degrees; 0 disables),
%   'preserve_boundary' (true), 'interpolate_point_data' (true) and
%   'record_parent_ids' (false). Roadmap 9.2; see doc/subdivide_surface.md.
%   Returns a struct: mesh, num_crease_edges, num_boundary_edges, warnings.

    o = mio.internal.opts(varargin, {'scheme', 'loop', 'levels', 1, 'crease_angle', 60, ...
        'preserve_boundary', true, 'interpolate_point_data', true, ...
        'record_parent_ids', false}, 'mio.subdivide_surface');
    r = meshioplusplus_mex('mio_subdivide_surface', mesh, mio.internal.str(o.scheme, 'scheme'), ...
        double(o.levels), double(o.crease_angle), logical(o.preserve_boundary), ...
        logical(o.interpolate_point_data), logical(o.record_parent_ids));
    if isempty(r.warnings)
        r.warnings = {};
    else
        r.warnings = strsplit(r.warnings, char(10));
    end
end
