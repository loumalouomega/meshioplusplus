function m = crop_predicate(mesh, array, varargin)
%CROP_PREDICATE  Keep the cells whose cell data ARRAY passes a comparison.
%   M = mio.crop_predicate(MESH, ARRAY, ...) keeps every cell with
%   ARRAY <compare> VALUE. Options: 'compare' ('<', '<=', '>', '>=', '==' or
%   the not-equal operator, '<' by default), 'value' (0) and 'record_ids'
%   (false). There is deliberately no 'mode': a cell-data predicate is
%   already one value per cell.

    o = mio.internal.opts(varargin, {'compare', '<', 'value', 0, 'record_ids', false}, ...
        'mio.crop_predicate');
    % The not-equal operator is spelt with an exclamation mark, built here
    % from its two characters.
    ops = {'<', '<=', '>', '>=', '==', [char(33) '=']};
    m = meshioplusplus_mex('mio_crop_predicate', mesh, mio.internal.str(array, 'array'), ...
        mio.internal.choice(o.compare, ops, 'compare'), double(o.value), logical(o.record_ids));
end
