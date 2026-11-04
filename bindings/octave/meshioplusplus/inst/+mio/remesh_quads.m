function r = remesh_quads(mesh, varargin)
%REMESH_QUADS  Remesh a triangle surface into field-aligned quads.
%   R = mio.remesh_quads(MESH, ...) remeshes the triangle surface MESH into
%   quads aligned with a smooth cross field (roadmap 7.7.3) by an installed
%   QuadriFlow or Instant Meshes, run as a separate process. No remesher
%   installed fails naming the programs. Options:
%     'target_faces' the number of faces asked for; 0 takes half the triangles (0)
%     'backend'      'auto' (QuadriFlow, else Instant Meshes) or a name of
%                    mio.remesh_quads_backends ('auto')
%     'sharp'        align with and keep the sharp edges (false)
%     'boundary'     keep the boundary of an open surface (true)
%     'adaptive'     let the face size follow the curvature (QuadriFlow) (false)
%     'dominant'     allow triangles where quads fit badly (Instant Meshes) (false)
%     'seed'         the program's random seed (0)
%     'max_cells'    refuse a target or result of more faces (5e7)
%   R is a struct: mesh, num_quads, num_triangles, num_points,
%   min_scaled_jacobian, mean_scaled_jacobian, hausdorff_to_surface and
%   backend.

    o = mio.internal.opts(varargin, {'target_faces', 0, 'backend', 'auto', 'sharp', false, ...
        'boundary', true, 'adaptive', false, 'dominant', false, 'seed', 0, ...
        'max_cells', 5e7}, 'mio.remesh_quads');
    r = meshioplusplus_mex('mio_remesh_quads', mesh, double(o.target_faces), ...
        mio.internal.str(o.backend, 'backend'), logical(o.sharp), logical(o.boundary), ...
        logical(o.adaptive), logical(o.dominant), double(o.seed), double(o.max_cells));
end
