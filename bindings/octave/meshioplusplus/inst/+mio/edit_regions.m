function m = edit_regions(mesh, op, inputs, varargin)
%EDIT_REGIONS  Combine, rename, retag or delete named regions; returns a new mesh.
%   OP is 'union', 'intersection', 'difference', 'rename', 'retag' or
%   'delete'; INPUTS a cellstr of region names. Options: 'output' (''),
%   'kind' ([] for any), 'dim' ([]), 'tag' ([]), 'keep_inputs' (true).

    o = mio.internal.opts(varargin, {'output', '', 'kind', [], 'dim', [], 'tag', [], ...
        'keep_inputs', true}, 'mio.edit_regions');
    m = meshioplusplus_mex('mio_edit_regions', mesh, mio.internal.str(op, 'op'), ...
        mio.internal.cellstr_arg(inputs, 'inputs'), mio.internal.str(o.output, 'output'), ...
        mio.internal.str(o.kind, 'kind'), o.dim, o.tag, logical(o.keep_inputs));
end
