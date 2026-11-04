function m = extract_surface(mesh, varargin)
%EXTRACT_SURFACE  The boundary surface of a volume mesh, as a new mesh.
%   S = mio.extract_surface(MESH) returns the boundary facets of MESH's
%   volume cells (triangles and quads) as a new mesh; MESH is unchanged.
%   Options: 'record_parent_ids' (false) attaches the cell data
%   'surface:parent_cell' naming each facet's source cell.

    o = mio.internal.opts(varargin, {'record_parent_ids', false}, 'mio.extract_surface');
    m = meshioplusplus_mex('mio_extract_surface', mesh, logical(o.record_parent_ids));
end
