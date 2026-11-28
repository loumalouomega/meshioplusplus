function out = convert_coordinates_crs(mesh, target_crs)
%CONVERT_COORDINATES_CRS  Reproject the points to another CRS.
%   OUT = mio.convert_coordinates_crs(MESH, TARGET_CRS) reprojects from the CRS
%   MESH records to TARGET_CRS (roadmap 9.10.2). Needs a build with PROJ
%   (mio.has_crs_provider()).

    out = meshioplusplus_mex('mio_convert_coordinates_crs', mesh, target_crs);
end
