function r = exaggerate(mesh, varargin)
%EXAGGERATE  Exaggerate a deviation map for a figure.
%   R = mio.exaggerate(MESH, ...) returns a copy of the surface MESH with each
%   point displaced along its normal by 'scale' times its deviation, clamped
%   to +-'clamp' (0: the suggested colour range) and smoothed over
%   'smooth_fraction' of the bounding-box diagonal (roadmap 7.4.4). Options:
%     'field'           the signed deviation point data ('deviation:distance')
%     'scale'           displacement per unit deviation (10)
%     'clamp'           clamp the deviation first, 0 automatic (0)
%     'smooth_fraction' smoothing radius over the diagonal (0.01)
%   R is a struct: mesh, clamp, radius, num_moved and warnings.

    o = mio.internal.opts(varargin, {'field', 'deviation:distance', 'scale', 10, ...
        'clamp', 0, 'smooth_fraction', 0.01}, 'mio.exaggerate');
    r = meshioplusplus_mex('mio_exaggerate', mesh, mio.internal.str(o.field, 'field'), ...
        double(o.scale), double(o.clamp), double(o.smooth_fraction));
    if isempty(r.warnings)
        r.warnings = cell(0, 1);
    else
        w = strsplit(r.warnings, char(10));
        r.warnings = w(:);
    end
end
