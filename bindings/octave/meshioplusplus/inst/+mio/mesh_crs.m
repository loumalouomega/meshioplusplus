function crs = mesh_crs(mesh)
%MESH_CRS  The coordinate reference system a mesh records.
%   CRS = mio.mesh_crs(MESH) is field_data 'meshio:crs' as text, '' when none.

    crs = meshioplusplus_mex('mio_mesh_crs', mesh);
end
