function m = platonic(varargin)
%PLATONIC  A Platonic solid's surface.
%   M = mio.platonic(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'name'       'tetrahedron', 'cube', 'octahedron', 'dodecahedron' or
%                  'icosahedron' ('icosahedron')
%     'radius'     circumradius (1)
%     'center'     centre of the shape ([0 0 0])
%     'max_cells'  refuse a mesh larger than this, by name (20000000)
%   Returns a mio.MeshHandle. See doc/primitives.md.

    o = mio.internal.opts(varargin, {'name', 'icosahedron', 'radius', 1, 'center', [0 0 0], ...
        'max_cells', 20000000}, 'mio.platonic');
    solids = {'tetrahedron', 'cube', 'octahedron', 'dodecahedron', 'icosahedron'};
    faces = [4, 6, 8, 12, 20];
    hit = [];
    if ischar(o.name)
        hit = find(strcmp(o.name, solids), 1);
    end
    if isempty(hit)
        error('meshioplusplus:argument', 'unknown platonic solid ''%s''', ...
              mio.internal.str(o.name, 'name'));
    end
    m = mio.internal.od_generate('platonic', faces(hit), 'origin', o.center, ...
        'radius', o.radius, 'max_cells', o.max_cells);
end
