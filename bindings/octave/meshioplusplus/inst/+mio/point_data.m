function x = point_data(mesh, name)
%POINT_DATA  A copy of point data NAME in its stored class (num_rows x components...).

    x = meshioplusplus_mex('mio_point_data', mesh, mio.internal.str(name, 'name'));
end
