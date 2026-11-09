function r = fit_primitives(mesh, varargin)
%FIT_PRIMITIVES  Fit a plane, sphere, cylinder, cone or torus to each patch.
%   Options: 'patch_data' ('segment:patch'; '' fits the whole surface),
%   'tolerance' (1e-3), 'relative' (true: a fraction of each patch's
%   bounding-box diagonal), 'min_inlier_fraction' (0.9), 'types' (cellstr of
%   the primitives to try; all five by default) and 'max_iterations' (50).
%   Roadmap 7.1: 'region' fits a named cell region as one patch, 'seeds' (an
%   n x 3 matrix) grows and fits a feature from each point, crossing edges
%   below 'feature_angle' (30), peeling rim faces turning more than
%   'rim_angle' (15), written as the cell region 'output_region'
%   ('fit:feature'), and 'used_points' ('auto', 'all', '3sigma', '2sigma',
%   '1sigma') trims outliers.
%   Returns a struct: mesh (with cell data 'fit:type', 'fit:rms', 'fit:max'
%   and field data 'fit:patch_*'), num_patches, num_plane, num_sphere,
%   num_cylinder, num_cone, num_torus, num_unfitted, num_feature_cells and
%   warnings.

    o = mio.internal.opts(varargin, {'patch_data', 'segment:patch', 'tolerance', 1e-3, ...
        'relative', true, 'min_inlier_fraction', 0.9, ...
        'types', {'plane', 'sphere', 'cylinder', 'cone', 'torus'}, 'max_iterations', 50, ...
        'region', '', 'seeds', [], 'feature_angle', 30, 'rim_angle', 15, ...
        'output_region', 'fit:feature', 'used_points', 'auto'}, 'mio.fit_primitives');
    used = find(strcmp({'auto', 'all', '3sigma', '2sigma', '1sigma'}, o.used_points));
    if isempty(used)
        error('meshioplusplus:argument', 'meshio++: fit_primitives: unknown used_points %s', ...
            o.used_points);
    end
    seeds = [];
    if ~isempty(o.seeds)
        seeds = reshape(double(o.seeds)', 1, []);
    end
    r = meshioplusplus_mex('mio_fit_primitives', mesh, ...
        mio.internal.str(o.patch_data, 'patch_data'), double(o.tolerance), ...
        logical(o.relative), double(o.min_inlier_fraction), ...
        types_mask(o.types, 'fit_primitives'), double(o.max_iterations), ...
        mio.internal.str(o.region, 'region'), seeds, double(o.feature_angle), ...
        double(o.rim_angle), mio.internal.str(o.output_region, 'output_region'), ...
        double(used - 1));
    if isfield(r, 'warnings')
        if isempty(r.warnings)
            r.warnings = cell(0, 1);
        else
            w = strsplit(r.warnings, char(10));
            r.warnings = w(:);
        end
    end
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
