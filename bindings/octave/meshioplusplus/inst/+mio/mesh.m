function m = mesh()
%MESH  A new, empty mesh (a mio.MeshHandle).

    m = meshioplusplus_mex('mio_mesh_create');
end
