function r = distance_to_surface(mesh, surface, varargin)
%DISTANCE_TO_SURFACE  Attach the signed distance to SURFACE on MESH.
%   R = mio.distance_to_surface(MESH, SURFACE) returns a struct: mesh (a copy
%   of MESH carrying the sdf:distance data array), num_banded and quality
%   (the verdict for SURFACE; see mio.surface_watertight_check). Options:
%     'sign'              'pseudonormal' (default), 'winding-number' or
%                         'unsigned'
%     'location'          'corner'/'point' (point data, default) or
%                         'center'/'centre'/'cell' (cell data)
%     'band'              clamp distances beyond this; 0 is the full field (0)
%     'record_inside'     attach the sdf:inside array (false)
%     'watertight_check'  'off', 'warn' (default) or 'error'

    o = mio.internal.opts(varargin, {'sign', 'pseudonormal', 'location', 'corner', ...
        'band', 0, 'record_inside', false, 'watertight_check', 'warn'}, ...
        'mio.distance_to_surface');
    sign = mio.internal.od_code(o.sign, {'unsigned', 'pseudonormal', 'winding-number'}, 'sign');
    location = mio.internal.od_sdf_location(o.location);
    check = mio.internal.od_code(o.watertight_check, {'off', 'warn', 'error'}, 'watertight check');
    r = meshioplusplus_mex('mio_distance_to_surface', mesh, surface, sign, location, ...
        double(o.band), logical(o.record_inside), check);
end
