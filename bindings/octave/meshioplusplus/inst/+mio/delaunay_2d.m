function r = delaunay_2d(mesh)
%DELAUNAY_2D  Delaunay triangulation of MESH's points in the x, y plane.
%   Every point is kept under one counter-clockwise triangle block. Returns
%   a struct: mesh, num_triangles, num_duplicates and num_boundary.

    r = meshioplusplus_mex('mio_delaunay_2d', mesh);
end
