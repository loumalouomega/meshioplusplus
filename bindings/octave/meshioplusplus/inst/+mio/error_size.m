function r = error_size(mesh, target_error, varargin)
%ERROR_SIZE  A per-point size field from a per-cell error estimate.
%   R = mio.error_size(MESH, TARGET_ERROR, ...) turns the per-cell error
%   indicator 'error_array' (what mio.estimate_error writes) into a size field
%   by the adaptive procedure of Zienkiewicz and Zhu (1987; roadmap 6.10.2):
%   with N cells evaluated, a cell of diameter h whose indicator is eta asks
%   h * (target_error / sqrt(N) / eta)^(1 / order), clamped to
%   [min_size, max_size] (0: no bound); a point takes the least of its cells'
%   (NaN where none is evaluated). The sizes are written as point data
%   'output'. TARGET_ERROR is positive. Options:
%     'error_array' the one-component cell-data indicator ('error:zz')
%     'order'       the indicator's rate in the cell size (1)
%     'min_size'    lower bound, 0 for none (0)
%     'max_size'    upper bound, 0 for none (0)
%     'output'      the point-data array written ('size')
%     'overwrite'   replace an existing array of that name (false)
%   R is a struct: mesh, num_cells (evaluated), num_skipped, num_clamped,
%   min_size and max_size (the least and greatest size written).

    o = mio.internal.opts(varargin, {'error_array', 'error:zz', 'order', 1, 'min_size', 0, ...
        'max_size', 0, 'output', 'size', 'overwrite', false}, 'mio.error_size');
    r = meshioplusplus_mex('mio_error_size', mesh, double(target_error), ...
        mio.internal.str(o.error_array, 'error_array'), double(o.order), ...
        double(o.min_size), double(o.max_size), mio.internal.str(o.output, 'output'), ...
        logical(o.overwrite));
end
