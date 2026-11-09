function r = offset_surface(mesh, distance, varargin)
%OFFSET_SURFACE  Offset a triangle and quad surface.
%   R = mio.offset_surface(MESH, DISTANCE, ...) moves each point along its
%   offset direction so every face's plane moves by DISTANCE (positive along
%   the faces' normals); 'smooth_iterations' (0) smooths the directions
%   (roadmap 7.8.4). R is a struct: mesh, num_not_visible, num_clamped,
%   num_folded, num_rim_quads, max_move and warnings.

    o = mio.internal.opts(varargin, {'smooth_iterations', 0}, 'mio.offset_surface');
    r = meshioplusplus_mex('mio_offset_surface', mesh, double(distance), ...
        double(o.smooth_iterations));
    r = mio.internal.split_warnings(r);
end
