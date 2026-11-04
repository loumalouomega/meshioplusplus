function names = remesh_quads_backends()
%REMESH_QUADS_BACKENDS  The quad remeshers mio.remesh_quads takes.
%   Besides 'auto': every registered quad remesher, the installed
%   'quadriflow' and 'instant_meshes' programs (found on PATH, or through
%   MESHIOPLUSPLUS_QUADRIFLOW and MESHIOPLUSPLUS_INSTANT_MESHES), sorted, as
%   an n x 1 cellstr.

    names = meshioplusplus_mex('mio_remesh_quads_backends');
end
