function r = mesh_boolean(a, b, varargin)
%MESH_BOOLEAN  Watertight surface boolean of two closed triangle meshes.
%   R = mio.mesh_boolean(A, B, ...) with 'op' ('union', 'difference',
%   'intersection') and 'backend' ([] for the registered provider) runs the
%   registered boolean provider (Manifold; roadmap 9.1). Fails when none is
%   registered (mio.has_boolean_provider()). See doc/mesh_boolean.md.
%   Returns a struct: mesh, num_triangles.

    o = mio.internal.opts(varargin, {'op', 'union', 'backend', []}, 'mio.mesh_boolean');
    r = meshioplusplus_mex('mio_mesh_boolean', a, b, mio.internal.str(o.op, 'op'), ...
        mio.internal.str(o.backend, 'backend'));
end
