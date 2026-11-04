function m = extract_skin(mesh, varargin)
%EXTRACT_SKIN  The outer skin of a mesh, as a new mesh.
%   S = mio.extract_skin(MESH) returns the skin facets of MESH.
%   Options: 'linearize' (false) emits corner nodes only (triangle/quad
%   output).

    o = mio.internal.opts(varargin, {'linearize', false}, 'mio.extract_skin');
    m = meshioplusplus_mex('mio_extract_skin', mesh, logical(o.linearize));
end
