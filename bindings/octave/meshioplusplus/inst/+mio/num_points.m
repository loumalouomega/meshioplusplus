function out = num_points(mesh)
%NUM_POINTS  Number of points.

    out = meshioplusplus_mex('mio_num_points', mesh);
end
