function r = hollow(mesh, wall, varargin)
%HOLLOW  Give a closed surface an inner wall.
%   R = mio.hollow(MESH, WALL, ...) adds the surface's offset by -WALL,
%   reversed, inside it; an open surface is refused (roadmap 7.8.4). R is as
%   mio.offset_surface's.

    o = mio.internal.opts(varargin, {'smooth_iterations', 0}, 'mio.hollow');
    r = meshioplusplus_mex('mio_hollow', mesh, double(wall), double(o.smooth_iterations));
    r = mio.internal.split_warnings(r);
end
