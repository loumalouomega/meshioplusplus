function names = tetrahedralize_backends()
%TETRAHEDRALIZE_BACKENDS  The backend names mio.tetrahedralize takes.
%   Besides 'auto': 'octree', then every registered tetrahedral mesher
%   ('netgen' when built in, the installed 'tetgen' and 'gmsh' programs,
%   found on PATH or through MESHIOPLUSPLUS_TETGEN and MESHIOPLUSPLUS_GMSH),
%   sorted, as an n x 1 cellstr.

    names = meshioplusplus_mex('mio_tetrahedralize_backends');
end
