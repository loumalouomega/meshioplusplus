function r = select(mesh, name, varargin)
%SELECT  Write a named region chosen by geometry, data, quality, component or growth.
%   R = mio.select(MESH, NAME, ...) selects points, cells or sides and writes
%   them as one named region; the mesh comes back unchanged plus that region
%   (a region of the same name and kind is replaced). A candidate is selected
%   when it passes every criterion given; at least one is required. See
%   doc/select.md. Options (defaults in parentheses):
%     'kind'              'cell', 'point' or 'side' (faces of the 3-D cells,
%                         or edges of the 2-D cells when there are none)
%                         ('cell')
%     'box'               [xmin ymin zmin xmax ymax zmax], or [] ([])
%     'plane'             [px py pz nx ny nz]: inside is
%                         (p - point) . normal >= 0 ([])
%     'sphere'            [cx cy cz r], the closed ball ([])
%     'where_array'       a scalar point data (points) or cell data (cells; a
%                         side takes its cell's value) array; '' is off ('')
%     'where_op'          '<', '<=', '>', '>=', '==' or '~=' ('<'; also lt, le, gt, ge, eq, ne)
%     'where_value'       the value compared against (0)
%     'quality_metric'    the same comparison on a quality metric (as
%                         mio.attach_quality names them, without the
%                         quality: prefix); cells and sides only ('')
%     'quality_op'        ('<')
%     'quality_value'     (0)
%     'component'         the k-th node-connected component, 1-based, in the
%                         order mio.split(..., 'by', 'component') returns its
%                         pieces; [] is off ([])
%     'seed'              [x y z]: grow a patch from the facet nearest it
%                         across edges no sharper than feature_angle; [] is
%                         off ([])
%     'feature_angle'     the largest dihedral angle, degrees, growth crosses
%                         (30)
%     'normal_tolerance'  the largest angle, degrees, between a grown facet's
%                         normal and the seed's; 0 is off (0)
%     'mode'              how a cell or side reduces its nodes for box, plane
%                         and sphere: 'all', 'any' or 'centroid' ('all')
%     'boundary_only'     side regions consider only facets owned by one cell
%                         (true)
%     'tag'               the region's integer tag; -1 for none (-1)
%   R is a struct: mesh, num_candidates, num_selected, num_grown (-1 without
%   a seed) and num_components (-1 without that criterion).

    o = mio.internal.opts(varargin, {'kind', 'cell', 'box', [], 'plane', [], 'sphere', [], ...
        'where_array', '', 'where_op', '<', 'where_value', 0, 'quality_metric', '', ...
        'quality_op', '<', 'quality_value', 0, 'component', [], 'seed', [], ...
        'feature_angle', 30, 'normal_tolerance', 0, 'mode', 'all', 'boundary_only', true, ...
        'tag', -1}, 'mio.select');
    r = meshioplusplus_mex('mio_select', mesh, mio.internal.str(name, 'name'), ...
        mio.internal.str(o.kind, 'kind'), double(o.tag), mio.internal.str(o.mode, 'mode'), ...
        logical(o.boundary_only), vec(o.box), vec(o.plane), vec(o.sphere), ...
        mio.internal.str(o.where_array, 'where_array'), op(o.where_op, 'where_op'), ...
        double(o.where_value), mio.internal.str(o.quality_metric, 'quality_metric'), ...
        op(o.quality_op, 'quality_op'), double(o.quality_value), vec(o.component), ...
        vec(o.seed), double(o.feature_angle), double(o.normal_tolerance));
end

function v = vec(x)
    % R's as.numeric(unlist(x)); [] stays [] (the criterion is off).
    if iscell(x)
        x = cell2mat(x);
    end
    if isempty(x)
        v = [];
    else
        v = double(x(:));
    end
end

function s = op(x, what)
    % MATLAB spells not-equal ~=; the C comparison table knows it as ne.
    s = mio.internal.str(x, what);
    if strcmp(s, '~=')
        s = 'ne';
    end
end
