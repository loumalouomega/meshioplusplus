function r = refine_patches(mesh, varargin)
%REFINE_PATCHES  Refine a segmentation until each patch is one surface.
%   Variational shape approximation on the 'segment:patch' labels from
%   mio.segment. Options: 'tolerance' (1e-3), 'relative' (true),
%   'normal_weight' (1), 'epsilon' (1e-3), 'feature_angle' (30),
%   'max_iterations' (10), 'min_patch_cells' (3), 'ransac_seed' ([] takes the
%   default; a whole number, exact below 2^53), 'types' ([] for every
%   primitive, else a cellstr), 'allow_split' (true), 'allow_merge' (true)
%   and 'prefix' ('patch'). Returns a struct: mesh, num_patches,
%   num_iterations, num_splits, num_merges, num_absorbed, num_primitive,
%   num_freeform, num_facets, tolerance, explained_area_fraction,
%   initial_explained_area_fraction, hausdorff_max, hausdorff_rms,
%   initial_energy, energy and stop ('unchanged', 'small_gain',
%   'iteration_cap' or 'energy_rose').

    o = mio.internal.opts(varargin, {'tolerance', 1e-3, 'relative', true, 'normal_weight', 1, ...
        'epsilon', 1e-3, 'feature_angle', 30, 'max_iterations', 10, 'min_patch_cells', 3, ...
        'ransac_seed', [], 'types', [], 'allow_split', true, 'allow_merge', true, ...
        'prefix', 'patch'}, 'mio.refine_patches');
    if isempty(o.types) && isnumeric(o.types)
        mask = 62;
    else
        mask = types_mask(o.types, 'refine_patches');
    end
    if isempty(o.ransac_seed)
        seed = [];
    else
        seed = double(o.ransac_seed);
    end
    r = meshioplusplus_mex('mio_refine_patches', mesh, double(o.tolerance), ...
        logical(o.relative), double(o.normal_weight), double(o.epsilon), ...
        double(o.feature_angle), double(o.max_iterations), double(o.min_patch_cells), seed, ...
        mask, logical(o.allow_split), logical(o.allow_merge), mio.internal.str(o.prefix, 'prefix'));
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
