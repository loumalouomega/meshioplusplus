function labels = partition_labels(mesh, nparts, varargin)
%PARTITION_LABELS  The part id of every cell, without building the pieces.
%   LABELS = mio.partition_labels(MESH, NPARTS, ...) is a num_cells x 1
%   column of part ids: ids, not indices, so deliberately NOT shifted to
%   1-based. Options as mio.partition: 'method' ('auto'), 'imbalance'
%   (0.03), 'partition_mode' ('eco'), 'seed' (0), 'weights' ([]).

    o = mio.internal.opts(varargin, {'method', 'auto', 'imbalance', 0.03, ...
        'partition_mode', 'eco', 'seed', 0, 'weights', []}, 'mio.partition_labels');
    labels = meshioplusplus_mex('mio_partition_labels', mesh, double(nparts), ...
        mio.internal.str(o.method, 'method'), double(o.imbalance), ...
        mio.internal.str(o.partition_mode, 'partition_mode'), double(o.seed), ...
        mio.internal.str(o.weights, 'weights'), mio.num_cells(mesh));
end
