function m = reference_cells(varargin)
%REFERENCE_CELLS  One reference cell of each requested type.
%   M = mio.reference_cells() builds every supported cell type, one block
%   each. Options:
%     'types'  a cellstr of meshio cell type names, or one comma-separated
%              string; '' (default) means every supported type
%     'pairs'  two cells sharing a facet instead of one (linear types; false)

    o = mio.internal.opts(varargin, {'types', '', 'pairs', false}, 'mio.reference_cells');
    types = o.types;
    if iscell(types)
        types = strjoin(mio.internal.cellstr_arg(types, 'types'), ',');
    end
    m = meshioplusplus_mex('mio_reference_cells', mio.internal.str(types, 'types'), ...
        logical(o.pairs));
end
