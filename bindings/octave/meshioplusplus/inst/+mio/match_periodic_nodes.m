function p = match_periodic_nodes(mesh, slave, master, varargin)
%MATCH_PERIODIC_NODES  Pair the nodes of region SLAVE with those of region MASTER.
%   Options: 'matrix' (4x4 affine mapping slave to master; [] for identity),
%   'translate' (3-vector, used when 'matrix' is []), 'atol' (1e-8),
%   'require_complete' (true). Returns a struct: slave, master, unmatched
%   (1-based columns), num_fixed and max_residual.

    o = mio.internal.opts(varargin, {'matrix', [], 'translate', [], 'atol', 1e-8, ...
        'require_complete', true}, 'mio.match_periodic_nodes');
    matrix = o.matrix;
    if isempty(matrix)
        matrix = eye(4);
        if ~isempty(o.translate)
            matrix(1:3, 4) = double(o.translate(:));
        end
    end
    p = meshioplusplus_mex('mio_match_periodic_nodes', mesh, mio.internal.str(slave, 'slave'), ...
        mio.internal.str(master, 'master'), double(matrix), double(o.atol), ...
        logical(o.require_complete));
end
