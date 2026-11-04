function m = region_adjacency(mesh, varargin)
%REGION_ADJACENCY  The conforming facets shared by Cell regions, as a mesh.
%   M = mio.region_adjacency(MESH) uses every Cell region (or, with fewer
%   than two, the cell blocks as groups); 'regions' (a cellstr, [] for all)
%   selects them. The mesh carries interface parent-cell/facet ids and
%   measure. See doc/region_adjacency.md.

    o = mio.internal.opts(varargin, {'regions', []}, 'mio.region_adjacency');
    m = meshioplusplus_mex('mio_region_adjacency', mesh, ...
        mio.internal.cellstr_arg(o.regions, 'regions'));
end
