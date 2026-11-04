function parts = partition(mesh, nparts, varargin)
%PARTITION  Partition a mesh into NPARTS pieces.
%   PARTS = mio.partition(MESH, NPARTS, ...). Options:
%     'method'         'sfc', 'kahip' or 'auto' ('auto')
%     'imbalance'      KaHIP only: allowed imbalance fraction in (0, 1) (0.03)
%     'partition_mode' KaHIP only: 'fast', 'eco' or 'strong' ('eco')
%     'seed'           KaHIP only: random seed (0)
%     'record_ids'     attach the original ids (false)
%     'ghost_layers'   grow each piece by this many shared-node BFS layers
%                      of other parts' cells, tagged 'partition:ghost' (0)
%     'weights'        name of a scalar cell-data array of per-cell weights ([])
%   PARTS is an n x 1 struct array with fields part_id (0-based, an id, not
%   an index), mesh, point_map (1-based) and cell_maps (one 1-based map per
%   block). Every piece owns its handle.

    o = mio.internal.opts(varargin, {'method', 'auto', 'imbalance', 0.03, ...
        'partition_mode', 'eco', 'seed', 0, 'record_ids', false, 'ghost_layers', 0, ...
        'weights', []}, 'mio.partition');
    parts = meshioplusplus_mex('mio_partition', mesh, double(nparts), ...
        mio.internal.str(o.method, 'method'), double(o.imbalance), ...
        mio.internal.str(o.partition_mode, 'partition_mode'), double(o.seed), ...
        logical(o.record_ids), double(o.ghost_layers), mio.internal.str(o.weights, 'weights'));
end
