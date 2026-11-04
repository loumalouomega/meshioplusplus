function r = estimate_normals(mesh, varargin)
%ESTIMATE_NORMALS  Normals of a point cloud by local plane fits.
%   Writes 'normals' point data (NaN where no plane fits). Options: 'k'
%   (16 neighbours), 'orientation' ('consistent', 'viewpoint' or 'none'),
%   'viewpoint' ([0 0 0]) and 'record_variation' (false). Returns a struct:
%   mesh, num_degenerate and num_components.

    o = mio.internal.opts(varargin, {'k', 16, 'orientation', 'consistent', ...
        'viewpoint', [0 0 0], 'record_variation', false}, 'mio.estimate_normals');
    r = meshioplusplus_mex('mio_estimate_normals', mesh, double(o.k), ...
        mio.internal.str(o.orientation, 'orientation'), double(o.viewpoint), ...
        logical(o.record_variation));
end
