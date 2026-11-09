function r = measure_element_primitive(mesh, name)
%MEASURE_ELEMENT_PRIMITIVE  A measured element as a fitted primitive.
%   R = mio.measure_element_primitive(MESH, NAME) is the element NAME as a
%   primitive type (1 plane, 2 sphere, 3 cylinder, 4 cone) and its ten
%   parameters, a deviation reference (roadmap 7.4.2). R is a struct: type
%   and params.

    r = meshioplusplus_mex('mio_measure_element_primitive', mesh, mio.internal.str(name, 'name'));
end
