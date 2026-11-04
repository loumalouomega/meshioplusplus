function r = fit_primitives(mesh, varargin)
%FIT_PRIMITIVES  Fit a plane, sphere, cylinder, cone or torus to each patch.
%   Options: 'patch_data' ('segment:patch'; '' fits the whole surface),
%   'tolerance' (1e-3), 'relative' (true: a fraction of each patch's
%   bounding-box diagonal), 'min_inlier_fraction' (0.9), 'types' (cellstr of
%   the primitives to try; all five by default) and 'max_iterations' (50).
%   Returns a struct: mesh (with cell data 'fit:type', 'fit:rms', 'fit:max'
%   and field data 'fit:patch_*'), num_patches, num_plane, num_sphere,
%   num_cylinder, num_cone, num_torus and num_unfitted.

    o = mio.internal.opts(varargin, {'patch_data', 'segment:patch', 'tolerance', 1e-3, ...
        'relative', true, 'min_inlier_fraction', 0.9, ...
        'types', {'plane', 'sphere', 'cylinder', 'cone', 'torus'}, 'max_iterations', 50}, ...
        'mio.fit_primitives');
    r = meshioplusplus_mex('mio_fit_primitives', mesh, ...
        mio.internal.str(o.patch_data, 'patch_data'), double(o.tolerance), ...
        logical(o.relative), double(o.min_inlier_fraction), ...
        types_mask(o.types, 'fit_primitives'), double(o.max_iterations));
end

function mask = types_mask(types, fname)
% The primitive bit mask: plane 2, sphere 4, cylinder 8, cone 16, torus 32.
    known = {'plane', 'sphere', 'cylinder', 'cone', 'torus'};
    if isempty(types)
        mask = 0;
        return;
    end
    types = mio.internal.cellstr_arg(types, 'types');
    bad = types(~ismember(types, known));
    if ~isempty(bad)
        error('meshioplusplus:argument', 'meshio++: %s: unknown primitive type %s', fname, bad{1});
    end
    [~, loc] = ismember(unique(types), known);
    mask = sum(2 .^ loc);
end
