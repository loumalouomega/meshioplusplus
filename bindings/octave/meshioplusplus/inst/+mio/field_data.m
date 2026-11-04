function x = field_data(mesh, name)
%FIELD_DATA  A copy of field data NAME in its stored class (num_rows x components...).

    x = meshioplusplus_mex('mio_field_data', mesh, mio.internal.str(name, 'name'));
end
