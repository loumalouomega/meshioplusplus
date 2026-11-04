function r = refine(mesh, varargin)
%REFINE  Refine a mesh uniformly or selectively.
%   R = mio.refine(MESH, ...). Options:
%     'levels'            how many times to apply the refinement templates (1)
%     'record_parent_ids' attach the parent-cell array (false)
%     'cells'             global (block-major) 1-based cell indices to refine ([])
%     'region'            name of a region to refine: a cell region selects
%                         its cells, a point region every cell with a node in
%                         it; a side region is an error ([])
%     'where_array'       a scalar cell-data array to threshold ([])
%     'where_op'          '<', '<=', '>', '>=', '==' or '!=' ('<')
%     'where_value'       the threshold; a non-finite value never matches (0)
%     'closure'           'redgreen' (local), 'propagate' (converges to
%                         uniform refinement of the component) or 'balanced'
%                         (keeps hanging nodes, 2:1 balance; the output is
%                         then not conforming) ('redgreen')
%     'record_levels'     attach the 'refine:level' cell data (false)
%     'record_hierarchy'  attach 'refine:cell_id'/'refine:parent_id', the
%                         persistent parent/child hierarchy across meshes;
%                         also forces 'refine:entity' (false)
%   At most one of 'cells', 'region' and 'where_array' may be given; with
%   none, every cell is refined. R is a struct: mesh, point_map (n x 1,
%   1-based) and cell_maps (one 1-based first-child map per input block).

    o = mio.internal.opts(varargin, {'levels', 1, 'record_parent_ids', false, 'cells', [], ...
        'region', [], 'where_array', [], 'where_op', '<', 'where_value', 0, ...
        'closure', 'redgreen', 'record_levels', false, 'record_hierarchy', false}, ...
        'mio.refine');
    r = meshioplusplus_mex('mio_refine', mesh, double(o.levels), logical(o.record_parent_ids), ...
        double(o.cells), mio.internal.str(o.region, 'region'), ...
        mio.internal.str(o.where_array, 'where_array'), mio.internal.str(o.where_op, 'where_op'), ...
        double(o.where_value), mio.internal.str(o.closure, 'closure'), ...
        logical(o.record_levels), logical(o.record_hierarchy));
end
