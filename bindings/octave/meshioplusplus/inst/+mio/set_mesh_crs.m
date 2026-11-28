function mesh = set_mesh_crs(mesh, crs)
%SET_MESH_CRS  Record a coordinate reference system on a mesh.
%   MESH = mio.set_mesh_crs(MESH, CRS) records CRS ('EPSG:<code>' or WKT2) as
%   field_data 'meshio:crs'; '' clears it. MESH is modified in place.

    meshioplusplus_mex('mio_set_mesh_crs', mesh, crs);
end
