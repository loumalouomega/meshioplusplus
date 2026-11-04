function mesh = add_field_data(mesh, name, data)
%ADD_FIELD_DATA  Add field data NAME (any shape; a vector is 1-D).

    meshioplusplus_mex('mio_add_field_data', mesh, mio.internal.str(name, 'name'), data);
end
