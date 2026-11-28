function u = field_unit(mesh, array)
%FIELD_UNIT  The declared unit of a data array.
%   U = mio.field_unit(MESH, ARRAY) is [scale_to_SI m kg s K A mol cd], or []
%   when ARRAY has no declared unit.

    u = meshioplusplus_mex('mio_field_unit', mesh, array);
end
