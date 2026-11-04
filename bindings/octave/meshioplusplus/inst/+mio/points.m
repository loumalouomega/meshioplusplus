function out = points(mesh)
%POINTS  The points, a num_points x dim double matrix (one point per row).

    out = meshioplusplus_mex('mio_points', mesh);
end
