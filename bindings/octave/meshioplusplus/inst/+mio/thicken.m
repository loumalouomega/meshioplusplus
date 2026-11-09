function r = thicken(mesh, thickness, varargin)
%THICKEN  Turn a surface into a closed shell.
%   R = mio.thicken(MESH, THICKNESS, ...) offsets the surface to the two
%   sides of a wall THICKNESS thick ('side': 'inward' (default), 'outward' or
%   'both'), reverses the lower layer and joins the boundary edges by quads;
%   cell data thicken:side (roadmap 7.8.4). R is as mio.offset_surface's.

    o = mio.internal.opts(varargin, {'side', 'inward', 'smooth_iterations', 0}, 'mio.thicken');
    s = mio.internal.choice(o.side, {'inward', 'outward', 'both'}, 'side');
    r = meshioplusplus_mex('mio_thicken', mesh, double(thickness), double(s), ...
        double(o.smooth_iterations));
    r = mio.internal.split_warnings(r);
end
