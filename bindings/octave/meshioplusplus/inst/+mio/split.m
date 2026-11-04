function [pieces, keys] = split(mesh, varargin)
%SPLIT  Split a mesh into pieces.
%   [PIECES, KEYS] = mio.split(MESH, ...) returns PIECES, an n x 1 cell array
%   of mio.MeshHandle objects, and KEYS, the n x 1 cellstr of their keys
%   (R's names). Every piece owns its handle. Options:
%     'by'       'type', 'component', or 'region'/'tag' ('type')
%     'tag_name' for by = 'region', the integer cell-data array name ([])

    o = mio.internal.opts(varargin, {'by', 'type', 'tag_name', []}, 'mio.split');
    s = meshioplusplus_mex('mio_split', mesh, mio.internal.str(o.by, 'by'), ...
        mio.internal.str(o.tag_name, 'tag_name'));
    pieces = s.meshes;
    keys = s.keys;
end
