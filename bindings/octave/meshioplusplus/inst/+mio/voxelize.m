function r = voxelize(mesh, varargin)
%VOXELIZE  A hexahedral voxel grid around a surface.
%   R = mio.voxelize(MESH, ...) voxelizes the surface MESH. Give exactly one
%   of 'resolution' and 'cell_size'. Options (defaults in parentheses):
%     'resolution'        cell counts [nx ny nz] ([]); a zero third entry
%                         voxelizes a 2-D skin in 'plane'
%     'cell_size'         cubic cell size; 0 means unset (0)
%     'bounds'            explicit [xlo ylo zlo xhi yhi zhi] ([])
%     'padding'           grow the box by this on every side (0)
%     'padding_relative'  grow the box by this fraction of its diagonal (0)
%     'fill'              'all', 'surface' or 'inside' ('all')
%     'sign'              'pseudonormal', 'winding-number' or 'unsigned'
%                         ('pseudonormal')
%     'attach_occupancy'  attach the voxel:occupancy array (false)
%     'max_cells'         refuse a grid larger than this, by name (20000000)
%     'watertight_check'  'off', 'warn' or 'error' ('warn')
%     'plane'             'auto', 'xy', 'xz' or 'yz': the plane of a 2-D skin
%                         ('auto')
%   R is a struct: mesh, dims, origin and spacing (1 x 3 rows) and
%   num_occupied.

    o = mio.internal.opts(varargin, {'resolution', [], 'cell_size', 0, 'bounds', [], ...
        'padding', 0, 'padding_relative', 0, 'fill', 'all', 'sign', 'pseudonormal', ...
        'attach_occupancy', false, 'max_cells', 20000000, 'watertight_check', 'warn', ...
        'plane', 'auto'}, 'mio.voxelize');
    fill = mio.internal.od_code(o.fill, {'all', 'surface', 'inside'}, 'fill');
    plane = mio.internal.od_code(o.plane, {'auto', 'xy', 'xz', 'yz'}, 'plane');
    sign = mio.internal.od_code(o.sign, {'unsigned', 'pseudonormal', 'winding-number'}, 'sign');
    check = mio.internal.od_code(o.watertight_check, {'off', 'warn', 'error'}, 'watertight check');
    r = meshioplusplus_mex('mio_voxelize', mesh, double(o.resolution), double(o.cell_size), ...
        double(o.bounds), double(o.padding), double(o.padding_relative), fill, sign, ...
        logical(o.attach_occupancy), double(o.max_cells), check, plane);
end
