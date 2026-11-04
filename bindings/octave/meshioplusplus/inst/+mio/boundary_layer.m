function r = boundary_layer(mesh, first_height, varargin)
%BOUNDARY_LAYER  Grow prismatic boundary layers from a surface.
%   R = mio.boundary_layer(MESH, FIRST_HEIGHT, ...) grows layers from the
%   triangle and quad surface MESH (roadmap 6.9): wedges from triangles and
%   hexahedra from quads along each point's most normal normal, the first
%   FIRST_HEIGHT thick and each next 'growth' times the one before. The wall
%   and the cap come out as faces, with 'boundary_layer:layer' cell data.
%   Options:
%     'num_layers'               how many layers (3)
%     'growth'                   each layer over the one before (1.2)
%     'reverse'                  grow against the faces' normals (false)
%     'collision_fraction'       of the distance a point's ray hits (0.45)
%     'smooth_iterations'        smoothing passes (5)
%     'max_shrink'               halvings before a point collapses (6)
%     'collapse'                 'points' or 'layers' ('points')
%     'normal_smooth_iterations' Laplacian passes over the growth
%                                directions, 0 to 1000 (0)
%     'front_check'              test the last layer's faces (true)
%   R is a struct: mesh, num_cells, num_wedges, num_hexahedra, num_pyramids,
%   num_tetrahedra, num_shrunk, num_collapsed, num_not_visible,
%   min_scaled_jacobian, min_thickness, max_thickness, cap_closed,
%   num_fewer_layers and num_front_crossings.

    o = mio.internal.opts(varargin, {'num_layers', 3, 'growth', 1.2, 'reverse', false, ...
        'collision_fraction', 0.45, 'smooth_iterations', 5, 'max_shrink', 6, ...
        'collapse', 'points', 'normal_smooth_iterations', 0, 'front_check', true}, ...
        'mio.boundary_layer');
    r = meshioplusplus_mex('mio_boundary_layer', mesh, double(first_height), ...
        double(o.num_layers), double(o.growth), logical(o.reverse), ...
        double(o.collision_fraction), double(o.smooth_iterations), double(o.max_shrink), ...
        mio.internal.str(o.collapse, 'collapse'), double(o.normal_smooth_iterations), ...
        logical(o.front_check));
end
