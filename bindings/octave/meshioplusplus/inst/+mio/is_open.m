function tf = is_open(mesh)
%IS_OPEN  True while MESH has not been released.

    tf = meshioplusplus_mex('mio_mesh_is_open', mesh);
end
