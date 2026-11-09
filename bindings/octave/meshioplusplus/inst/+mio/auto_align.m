function r = auto_align(mesh, varargin)
%AUTO_ALIGN  Pose a scan from itself alone.
%   On its largest flat face, or a cylindrical part on the end of its axis,
%   from an area-weighted histogram of its facet normals (roadmap 7.3.5).
%   Options: 'angle_tolerance' (10 degrees), 'subdivisions' (3), 'region'
%   ('') and 'rotate_vector_data' (false). Returns a struct: mesh, matrix,
%   mode ('face' or 'axis'), primary, secondary, share, face_share,
%   axis_share and opening_area.

    o = mio.internal.opts(varargin, {'angle_tolerance', 10, 'subdivisions', 3, ...
        'region', '', 'rotate_vector_data', false}, 'mio.auto_align');
    r = meshioplusplus_mex('mio_auto_align', mesh, double(o.angle_tolerance), ...
        double(o.subdivisions), mio.internal.str(o.region, 'region'), ...
        logical(o.rotate_vector_data));
end
