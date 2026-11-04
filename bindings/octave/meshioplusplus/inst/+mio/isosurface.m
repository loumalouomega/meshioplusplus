function m = isosurface(mesh, array, isovalues, varargin)
%ISOSURFACE  Level sets of the point data ARRAY at ISOVALUES, as a new mesh.
%   The cells carry 'iso:value' and 'iso:index' cell data. Options:
%   'component' (-1: the row magnitude of a multi-component array, else a
%   0-based component) and 'record_parent_ids' (false).

    o = mio.internal.opts(varargin, {'component', -1, 'record_parent_ids', false}, ...
        'mio.isosurface');
    m = meshioplusplus_mex('mio_isosurface', mesh, mio.internal.str(array, 'array'), ...
        double(isovalues), double(o.component), logical(o.record_parent_ids));
end
