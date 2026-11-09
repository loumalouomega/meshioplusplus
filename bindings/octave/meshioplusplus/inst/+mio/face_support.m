function r = face_support(solid, scan, varargin)
%FACE_SUPPORT  The scan's support of each face of a fitted solid.
%   R = mio.face_support(SOLID, SCAN, ...) classes each face of a fit_cad
%   re-mesh as measured, partial or inferred by the share of its area within
%   'tolerance' (0: 1e-3 of the scan's diagonal) of the scan (roadmap
%   7.7.1). Options: 'measured_share' (0.9), 'inferred_share' (0.1) and
%   'hole_data' ('repair:hole'). R is a struct: mesh, num_measured,
%   num_partial, num_inferred, num_faces, share_within and tolerance.

    o = mio.internal.opts(varargin, {'tolerance', 0, 'measured_share', 0.9, ...
        'inferred_share', 0.1, 'hole_data', 'repair:hole'}, 'mio.face_support');
    r = meshioplusplus_mex('mio_face_support', solid, scan, double(o.tolerance), ...
        double(o.measured_share), double(o.inferred_share), ...
        mio.internal.str(o.hole_data, 'hole_data'));
end
