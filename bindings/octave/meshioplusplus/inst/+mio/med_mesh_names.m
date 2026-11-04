function names = med_mesh_names(path)
%MED_MESH_NAMES  The mesh names of a MED file, in file link order (cellstr).

    names = meshioplusplus_mex('mio_med_mesh_names', mio.internal.str(path, 'path'));
end
