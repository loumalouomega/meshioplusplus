function r = segment(mesh, varargin)
%SEGMENT  Split a surface into smooth patches by region growing.
%   Options: 'feature_angle' (30), 'normal_tolerance' (0: off),
%   'curvature_tolerance' (0: off), 'min_patch_cells' (1), 'region' (''),
%   'prefix' ('patch': regions are named '<prefix>:<k>'), 'method' ('angle'
%   or 'primitive') and 'tolerance' (1e-3, 'primitive' only). Returns a
%   struct: mesh (with cell data 'segment:patch', field data
%   'segment:adjacency' and one cell region per patch), num_patches,
%   num_merged, num_feature_edges and num_boundary_edges.

    o = mio.internal.opts(varargin, {'feature_angle', 30, 'normal_tolerance', 0, ...
        'curvature_tolerance', 0, 'min_patch_cells', 1, 'region', '', 'prefix', 'patch', ...
        'method', 'angle', 'tolerance', 1e-3}, 'mio.segment');
    r = meshioplusplus_mex('mio_segment', mesh, double(o.feature_angle), ...
        double(o.normal_tolerance), double(o.curvature_tolerance), double(o.min_patch_cells), ...
        mio.internal.str(o.region, 'region'), mio.internal.str(o.prefix, 'prefix'), ...
        mio.internal.str(o.method, 'method'), double(o.tolerance));
end
