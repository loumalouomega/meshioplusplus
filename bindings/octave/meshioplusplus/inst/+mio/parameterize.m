function r = parameterize(mesh, varargin)
%PARAMETERIZE  A UV atlas of a triangle surface.
%   Cuts the surface into charts, flattens and packs them into the unit
%   square; the result is split at the seams (point data 'uv' and
%   'parameterize:parent_point', cell data 'parameterize:chart'). Options:
%   'padding' (0.02), 'max_depth' (8) and 'max_stretch' (2). Returns a
%   struct: mesh, num_charts, num_flipped, stretch, packing_efficiency and
%   num_seam_points.

    o = mio.internal.opts(varargin, {'padding', 0.02, 'max_depth', 8, 'max_stretch', 2}, ...
        'mio.parameterize');
    r = meshioplusplus_mex('mio_parameterize', mesh, double(o.padding), double(o.max_depth), ...
        double(o.max_stretch));
end
