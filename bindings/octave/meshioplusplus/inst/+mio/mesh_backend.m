function s = mesh_backend()
%MESH_BACKEND  The native mesh backend: 'meshio', 'native' or 'kratos'.

    s = meshioplusplus_mex('mio_mesh_backend');
end
