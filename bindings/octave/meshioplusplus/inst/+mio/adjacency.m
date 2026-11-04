function g = adjacency(mesh, varargin)
%ADJACENCY  Adjacency graph in CSR form: a struct with indptr and indices (1-based).
%   Option: 'kind' ('node').

    o = mio.internal.opts(varargin, {'kind', 'node'}, 'mio.adjacency');
    g = meshioplusplus_mex('mio_adjacency', mesh, mio.internal.str(o.kind, 'kind'));
end
