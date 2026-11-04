function r = diff(a, b, varargin)
%DIFF  A structured comparison of two meshes.
%   R = mio.diff(A, B, ...). Options: 'atol' (0), 'rtol' (0), 'unordered'
%   (match points by spatial proximity instead of by index; false).
%   R is a struct:
%     verdict  'identical', 'equal_within_tolerance' or 'different'
%     points   struct: max_abs_error, max_rel_error, worst_index (1-based
%              flat index into the coordinates, 0 for none), num_exceeding
%     blocks   struct of n x 1 columns, one entry per block: type_mismatch,
%              count_mismatch (logical) and conn_mismatch_count

    o = mio.internal.opts(varargin, {'atol', 0, 'rtol', 0, 'unordered', false}, 'mio.diff');
    r = meshioplusplus_mex('mio_diff', a, b, double(o.atol), double(o.rtol), ...
        logical(o.unordered));
end
