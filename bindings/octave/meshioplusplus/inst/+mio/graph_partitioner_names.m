function names = graph_partitioner_names()
%GRAPH_PARTITIONER_NAMES  The registered graph partitioners.
%   NAMES = mio.graph_partitioner_names() is a cell array of the names
%   mio.partition 'method' accepts beyond 'sfc', 'kahip' and 'auto'.

    text = meshioplusplus_mex('mio_graph_partitioner_names');
    if isempty(text)
        names = {};
    else
        names = strsplit(text, char(10));
    end
end
