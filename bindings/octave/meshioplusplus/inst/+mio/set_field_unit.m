function mesh = set_field_unit(mesh, array, unit)
%SET_FIELD_UNIT  Declare the physical unit of a point or cell data array.
%   MESH = mio.set_field_unit(MESH, ARRAY, UNIT) stores UNIT ('m', 'mm',
%   'MPa', 'kg/m^3', ...) as field_data 'meshio:unit:<ARRAY>' (roadmap
%   9.10.1); MESH is modified in place. See doc/units.md.

    meshioplusplus_mex('mio_set_field_unit', mesh, array, unit);
end
