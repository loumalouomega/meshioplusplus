function med_write_multi(path, meshes, names, varargin)
%MED_WRITE_MULTI  Write several meshes into one MED file.
%   MESHES is a nonempty cell array of meshes, NAMES a cellstr with one unique
%   name per mesh. Option: 'version' ('4.1.0').

    o = mio.internal.opts(varargin, {'version', '4.1.0'}, 'mio.med_write_multi');
    if ~iscell(meshes)
        meshes = num2cell(meshes);
    end
    meshioplusplus_mex('mio_med_write_multi', mio.internal.str(path, 'path'), meshes, ...
        mio.internal.cellstr_arg(names, 'names'), mio.internal.str(o.version, 'version'));
end
