function release(mesh)
%RELEASE  Free MESH now. Idempotent; later calls on MESH raise.

    meshioplusplus_mex('mio_mesh_release', mesh);
end
