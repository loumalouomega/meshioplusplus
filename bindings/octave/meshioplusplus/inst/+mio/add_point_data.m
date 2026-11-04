function mesh = add_point_data(mesh, name, data)
%ADD_POINT_DATA  Add point data NAME: one row per point; the class (double, single,
%   int32, uint8, ...) is kept.

    meshioplusplus_mex('mio_add_point_data', mesh, mio.internal.str(name, 'name'), data);
end
