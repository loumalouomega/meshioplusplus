function r = clip(mesh, origin, normal, varargin)
%CLIP  Cut a surface by a plane, keep a side and cap it.
%   R = mio.clip(MESH, ORIGIN, NORMAL, ...) keeps the side below the plane
%   ('keep_below', true; false keeps the side above), splits the straddling
%   triangles and with 'cap' (true) triangulates the cut's closed loops
%   (roadmap 7.8.3). R is a struct: mesh, num_kept, num_split,
%   num_cap_triangles, num_loops, num_open_chains and warnings.

    o = mio.internal.opts(varargin, {'keep_below', true, 'cap', true}, 'mio.clip');
    r = meshioplusplus_mex('mio_clip', mesh, double(origin), double(normal), ...
        logical(o.keep_below), logical(o.cap));
    if isempty(r.warnings)
        r.warnings = cell(0, 1);
    else
        r.warnings = strsplit(r.warnings, char(10));
        r.warnings = r.warnings(:);
    end
end
