function r = octree_mesh(mesh, varargin)
%OCTREE_MESH  Mesh the inside of a skin on a balanced octree.
%   R = mio.octree_mesh(MESH, ...) builds conforming graded tetrahedra by
%   default ('output' 'tet'), or the castellated leaves ('castellated') with
%   surrogate-boundary regions and sdf:distance / sdf:closest_point; 'hex'
%   is the all-hexahedral dual. MESH is the skin (any surface; the boundary
%   of a volume mesh is taken). Every option left [] keeps the default
%   mio_octree_mesh_opts_init sets on the C side; those defaults are shown.
%   See doc/octree_mesh.md for every option and report field.
%     'output'               'tet' (default), 'castellated' or 'hex'
%     'min_size'             smallest leaf edge; 0 (default) is max_size / 32
%     'max_size'             largest leaf edge; 0 (default) is an eighth of
%                            the skin's bounding-box diagonal
%     'max_level'            cap on the finest level, at most 19; -1
%                            (default) from the sizes
%     'bounds'               [xlo ylo zlo xhi yhi zhi] to mesh instead of the
%                            padded bounding box
%     'padding'              padding on every side, world units (0)
%     'padding_relative'     padding as a fraction of the diagonal (0.1)
%     'curvature_angle'      degrees a leaf edge may subtend on the skin (15);
%                            0 disables the curvature size
%     'feature_angle'        edges sharper than this dihedral, in degrees,
%                            are features (30)
%     'corner_angle'         a feature turning by more than this is a corner
%                            (30)
%     'feature_size'         size along feature edges; 0 (default) is min_size
%     'gap_cells'            leaves across the local thickness (2); 0
%                            disables it
%     'gap_outside'          also resolve gaps outside the skin (false)
%     'size_array'           a point data array of target sizes on the skin
%     'surface_region'       restrict the skin to this cell region
%     'grading_layers'       rings of leaves around every cut leaf at its
%                            level (1)
%     'skin_grading'         'graded' (default) or 'uniform'
%     'surface_size'         the skin size under uniform grading; 0 (default)
%                            is the smallest the field asks for
%     'winding_threshold'    inside where the winding number exceeds this (0.5)
%     'winding_beta'         fast winding number accuracy (2); 0 sums every
%                            triangle
%     'conformity'           'stuffing' (default) or 'cleaving' (tet output)
%     'alpha_long'           stuffing threshold on long lattice edges (0.24799)
%     'alpha_short'          stuffing threshold on short edges (0.40173)
%     'preserve_features'    snap the boundary to feature corners and edges
%                            (true)
%     'fit_iterations'       rounds of optimisation (5); 0 disables them
%     'dihedral_floor_guard' undo optimisation moves that worsen the dihedral
%                            angles (true)
%     'fill'                 castellated output: 'inside' (default), 'all' or
%                            'cover'
%     'keep_outside'         tet output: also mesh the outside (false)
%     'attach_distance'      attach sdf:distance (true)
%     'attach_quality'       attach the quality:* cell data (false)
%     'pillow'               hex output: a boundary layer of hexahedra (true)
%     'hex_grading'          hex output: 'uniform' (default, one grid at the
%                            skin's size) or 'graded' (the octree's grading,
%                            through dual transition schemes)
%     'untangle_iterations'  hex output: untangling iterations (50)
%     'max_cells'            refuse an octree with more leaves (20000000)
%     'max_tets'             refuse an output with more cells (100000000)
%     'volume_fix'           tet output: run the volume fix after the
%                            optimisation (true)
%     'zero_volume_factor'   degenerate under this fraction of the longest
%                            edge cubed (1e-6)
%     'cap_flatness'         a flat tetrahedron is a cap under this height
%                            over its face, as a fraction of its size (0.01)
%     'crossing_margin'      a sliver's edges cross this far inside both, a
%                            fraction in [0, 0.5) (0.01)
%     'crossing_gap'         ... and come this close, relative to the
%                            shorter edge (0.05)
%   R is a struct: mesh, the report fields of mio_octree_mesh_report (output
%   as its name, guarantee_applies as a logical, origin and root_cells as
%   1 x 3 rows, every count as a double), leaves_per_level (an n x 1 column
%   from level 0), timings (an n x 2 cell of stage name and seconds) and
%   warnings, an n x 1 cellstr with one entry per warning.

    names = {'output', 'min_size', 'max_size', 'max_level', 'bounds', 'padding', ...
        'padding_relative', 'curvature_angle', 'feature_angle', 'corner_angle', ...
        'feature_size', 'gap_cells', 'gap_outside', 'size_array', 'surface_region', ...
        'grading_layers', 'skin_grading', 'surface_size', 'winding_threshold', ...
        'winding_beta', 'conformity', 'alpha_long', 'alpha_short', 'preserve_features', ...
        'fit_iterations', 'dihedral_floor_guard', 'fill', 'keep_outside', ...
        'attach_distance', 'attach_quality', 'pillow', 'hex_grading', ...
        'untangle_iterations', 'max_cells', 'max_tets', 'volume_fix', ...
        'zero_volume_factor', 'cap_flatness', 'crossing_margin', 'crossing_gap'};
    spec = cell(1, 2 * numel(names));
    spec(1:2:end) = names;
    spec(2:2:end) = {[]};
    o = mio.internal.opts(varargin, spec, 'mio.octree_mesh');
    r = meshioplusplus_mex('mio_octree_mesh', mesh, ...
        code(o.output, {'castellated', 'tet', 'hex'}, 'output'), ...
        code(o.conformity, {'stuffing', 'cleaving'}, 'conformity'), ...
        code(o.skin_grading, {'graded', 'uniform'}, 'skin grading'), ...
        code(o.fill, {'all', 'inside', 'cover'}, 'fill'), ...
        double(o.bounds), name(o.size_array, 'size_array'), ...
        name(o.surface_region, 'surface_region'), ...
        num(o.min_size), num(o.max_size), num(o.padding), num(o.padding_relative), ...
        num(o.curvature_angle), num(o.feature_angle), num(o.corner_angle), ...
        num(o.feature_size), num(o.gap_cells), num(o.surface_size), ...
        num(o.winding_threshold), num(o.winding_beta), num(o.alpha_long), ...
        num(o.alpha_short), num(o.max_level), num(o.grading_layers), ...
        num(o.fit_iterations), num(o.untangle_iterations), num(o.max_cells), ...
        num(o.max_tets), flag(o.gap_outside), flag(o.preserve_features), ...
        flag(o.dihedral_floor_guard), flag(o.keep_outside), flag(o.attach_distance), ...
        flag(o.attach_quality), flag(o.pillow), ...
        code(o.hex_grading, {'uniform', 'graded'}, 'hex grading'), ...
        flag(o.volume_fix), num(o.zero_volume_factor), num(o.cap_flatness), ...
        num(o.crossing_margin), num(o.crossing_gap));
    % One newline-separated string from the C side, one cellstr entry per
    % warning here.
    if isempty(r.warnings)
        r.warnings = cell(0, 1);
    else
        w = strsplit(r.warnings, sprintf('\n')).';
        if isempty(w{end})
            w(end) = []; % R's strsplit drops a trailing empty piece
        end
        r.warnings = w;
    end
    if isempty(r.stage_names)
        stages = cell(0, 1);
    else
        stages = strsplit(r.stage_names, sprintf('\n')).';
    end
    r.timings = [stages, num2cell(r.stage_seconds)];
    r = rmfield(r, {'stage_names', 'stage_seconds'});
end

function c = code(value, table, what)
    if isempty(value) && ~ischar(value)
        c = [];
        return;
    end
    c = mio.internal.od_code(value, table, what);
end

function x = num(x)
    if ~isempty(x)
        x = double(x);
    end
end

function x = flag(x)
    if ~isempty(x)
        x = logical(x);
    end
end

function x = name(x, what)
    x = mio.internal.str(x, what);
end
