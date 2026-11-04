function r = octree(mesh, varargin)
%OCTREE  A balanced octree (or quadtree) refined around a skin.
%   R = mio.octree(MESH, ...) builds the octree leaves around the skin MESH.
%   Options (defaults in parentheses):
%     'root_cells'           root cells per axis [nx ny nz] ([] is a single
%                            root cell)
%     'min_level'            every leaf is split down to at least this (0)
%     'max_level'            no leaf beyond this, at most 20 (4)
%     'bounds'               explicit [xlo ylo zlo xhi yhi zhi] ([])
%     'padding'              grow the box by this on every side (0)
%     'padding_relative'     grow the box by this fraction of its diagonal (0)
%     'dimension'            2 for a planar line skin (a quadtree), 3 for a
%                            surface, 0 to decide from the cell types (0)
%     'sample_distance'      refine further within this distance of the skin (0)
%     'balance'              'none', 'face' or 'full' 2:1 balance ('face')
%     'fill'                 'all', 'surface' or 'inside' ('all')
%     'sign'                 'pseudonormal', 'winding-number' or 'unsigned'
%                            ('pseudonormal')
%     'attach_distance'      attach the signed distance at leaf centres (false)
%     'even_odd'             even-odd rule instead of the winding number, 2-D
%                            (false)
%     'max_cells'            refuse more leaves than this (20000000)
%     'watertight_check'     'off', 'warn' or 'error' ('warn')
%     'transition'           'hanging', 'polyhedra' or 'templates': the output
%                            where a fine leaf meets a coarse one ('hanging')
%     'snap'                 move the outer nodes onto the skin without
%                            inverting a leaf (false)
%     'surface_level'        split cut leaves down to this level; negative is
%                            max_level (-1)
%     'curvature_threshold'  refinement criterion (0)
%     'feature_angle'        refinement criterion (0)
%     'size_array'           refinement criterion: a point data array of
%                            target sizes ('')
%     'plane'                'auto', 'xy', 'xz' or 'yz': the plane of a 2-D
%                            skin ('auto')
%   R is a struct: mesh, dimension, num_leaves, num_hanging, finest_level,
%   origin (1 x 3), root_size, num_snapped and num_refused.

    o = mio.internal.opts(varargin, {'root_cells', [], 'min_level', 0, 'max_level', 4, ...
        'bounds', [], 'padding', 0, 'padding_relative', 0, 'dimension', 0, ...
        'sample_distance', 0, 'balance', 'face', 'fill', 'all', 'sign', 'pseudonormal', ...
        'attach_distance', false, 'even_odd', false, 'max_cells', 20000000, ...
        'watertight_check', 'warn', 'transition', 'hanging', 'snap', false, ...
        'surface_level', -1, 'curvature_threshold', 0, 'feature_angle', 0, ...
        'size_array', '', 'plane', 'auto'}, 'mio.octree');
    transition = mio.internal.od_code(o.transition, {'hanging', 'polyhedra', 'templates'}, ...
        'transition');
    plane = mio.internal.od_code(o.plane, {'auto', 'xy', 'xz', 'yz'}, 'plane');
    balance = mio.internal.od_code(o.balance, {'none', 'face', 'full'}, 'balance');
    fill = mio.internal.od_code(o.fill, {'all', 'surface', 'inside'}, 'fill');
    sign = mio.internal.od_code(o.sign, {'unsigned', 'pseudonormal', 'winding-number'}, 'sign');
    check = mio.internal.od_code(o.watertight_check, {'off', 'warn', 'error'}, 'watertight check');
    r = meshioplusplus_mex('mio_octree', mesh, double(o.root_cells), double(o.min_level), ...
        double(o.max_level), double(o.bounds), double(o.padding), ...
        double(o.padding_relative), double(o.dimension), double(o.sample_distance), ...
        balance, fill, sign, logical(o.attach_distance), logical(o.even_odd), ...
        double(o.max_cells), check, transition, logical(o.snap), double(o.surface_level), ...
        double(o.curvature_threshold), double(o.feature_angle), ...
        mio.internal.str(o.size_array, 'size_array'), plane);
end
