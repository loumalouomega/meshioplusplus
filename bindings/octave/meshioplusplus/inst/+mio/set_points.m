function mesh = set_points(mesh, points)
%SET_POINTS  Set the points from a num_points x dim matrix (copied as double).

    meshioplusplus_mex('mio_set_points', mesh, points);
end
