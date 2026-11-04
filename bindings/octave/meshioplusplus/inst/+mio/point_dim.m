function out = point_dim(mesh)
%POINT_DIM  Point dimension (2 or 3).

    out = meshioplusplus_mex('mio_point_dim', mesh);
end
