function r = locate_points(mesh, points, varargin)
%LOCATE_POINTS  The cell of MESH containing each point.
%   R = mio.locate_points(MESH, POINTS) reduces MESH to simplices (as
%   mio.interpolate does) and finds, for each row of the n x 3 matrix POINTS,
%   the first containing simplex in ascending order; a point on a shared face
%   belongs to the lowest-numbered cell. Volume meshes and planar triangle
%   meshes are located. Option: 'tolerance' (1e-12): a barycentric weight at
%   least -tolerance counts as inside.
%   Returns a struct: cell (n x 1, the 1-based input cell, 0 outside), nodes
%   (n x 4, the 1-based point ids of the containing simplex, 0 for none or
%   padding), barycentric (n x 4) and num_located.

    o = mio.internal.opts(varargin, {'tolerance', 1e-12}, 'mio.locate_points');
    r = meshioplusplus_mex('mio_locate_points', mesh, double(points), double(o.tolerance));
end
