function r = profile_solid(sketch, features, varargin)
%PROFILE_SOLID  Extrude or revolve a sketch into a solid (a CAD provider).
%   R = mio.profile_solid(SKETCH, FEATURES, ...) builds a fit_curves sketch
%   into a solid by FEATURES, ';'-separated and applied in order:
%   'extrude:D[:symmetric]' (0 runs to the 'scan''s extent),
%   'revolve:OX,OY,OZ:DX,DY,DZ[:ANGLE]', 'fillet:R', 'chamfer:D', 'shell:W'
%   (roadmap 7.7.4). Options: 'step_path' (''), 'scan' ([]), 'provider'
%   ('') and 'deflection' (0). R is a struct: mesh, volume, area, num_faces,
%   valid and warnings.

    o = mio.internal.opts(varargin, {'step_path', '', 'scan', [], 'provider', '', ...
        'deflection', 0}, 'mio.profile_solid');
    r = meshioplusplus_mex('mio_profile_solid', sketch, mio.internal.str(features, 'features'), ...
        mio.internal.str(o.step_path, 'step_path'), o.scan, ...
        mio.internal.str(o.provider, 'provider'), double(o.deflection));
    r = mio.internal.split_warnings(r);
end
