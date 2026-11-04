function s = write_graph(mesh, path, varargin)
%WRITE_GRAPH  Write the adjacency graph to PATH. Options: 'kind' ('node'), 'format'
%   ([] to infer). Returns num_vertices, num_edges and max_degree.

    o = mio.internal.opts(varargin, {'kind', 'node', 'format', []}, 'mio.write_graph');
    s = meshioplusplus_mex('mio_write_graph', mesh, mio.internal.str(path, 'path'), ...
        mio.internal.str(o.kind, 'kind'), mio.internal.str(o.format, 'format'));
end
