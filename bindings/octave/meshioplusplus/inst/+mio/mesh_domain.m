function r = mesh_domain(boundary, varargin)
%MESH_DOMAIN  Triangulate a planar domain with MMG.
%   R = mio.mesh_domain(BOUNDARY, ...) triangulates the regions the 'line'
%   cells of BOUNDARY enclose by MMG2D_mmg2dmesh (roadmap 3.2.3,
%   doc/mmg.md), each a cell region 'subdomain:<k>'; 'holes' gives one point
%   inside each region to leave out. Options:
%     'hsiz'  a constant target edge size, negative for none (-1)
%     'size'  a one-component point-data name: the size per point ('')
%     'holes' [] or a k x 2 matrix of points (x, y), one per hole ([])
%     'hmin', 'hmax', 'hausd', 'hgrad', 'angle_detection', 'frozen',
%     'required' as in mio.adapt
%   R is the struct mio.adapt returns, plus num_subdomains and num_holes.

    o = mio.internal.opts(varargin, {'hsiz', -1, 'size', '', 'holes', [], 'hmin', -1, ...
        'hmax', -1, 'hausd', -1, 'hgrad', -1, 'angle_detection', 45, 'frozen', [], ...
        'required', {}}, 'mio.mesh_domain');
    holes = o.holes;
    if ~isempty(holes)
        holes = reshape(double(holes), [], 2);
    end
    r = meshioplusplus_mex('mio_mesh_domain', boundary, double(o.hsiz), ...
        mio.internal.str(o.size, 'size'), holes, double(o.hmin), double(o.hmax), ...
        double(o.hausd), double(o.hgrad), double(o.angle_detection), double(o.frozen), ...
        mio.internal.cellstr_arg(o.required, 'required'));
end
