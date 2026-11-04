function r = compute_sdf(mesh, varargin)
%COMPUTE_SDF  Generate a grid around a surface and its signed distance field.
%   R = mio.compute_sdf(MESH, ...) builds a voxel grid or an octree around
%   the surface MESH and samples the signed distance on it, in one call.
%   Options (defaults in parentheses):
%     'structure'         'voxel' or 'octree' ('voxel')
%     'resolution'        voxel cell counts [nx ny nz] ([]); give this or
%                         'cell_size' for a voxel grid -- an octree's finest
%                         cell is set by root_resolution and max_depth
%     'cell_size'         cubic voxel size; 0 means unset (0)
%     'bounds'            explicit [xlo ylo zlo xhi yhi zhi] ([])
%     'padding'           grow the box by this on every side (0)
%     'padding_relative'  grow the box by this fraction of its diagonal (0.1)
%     'root_resolution'   octree: cell count per axis of the root lattice (8)
%     'max_depth'         octree: how many refinement passes (4)
%     'band_cells'        octree: refine while |distance| <= this times the
%                         cell's own diagonal (1)
%     'record_levels'     octree: attach refine:level (true)
%     'max_cells'         refuse above this many cells; re-checked after
%                         every octree pass (20000000)
%     'sign'              'pseudonormal', 'winding-number' or 'unsigned'
%                         ('pseudonormal')
%     'location'          'corner' (or 'point') or 'center' (or 'centre',
%                         'cell'): where the field is sampled ('corner')
%     'band'              clamp distances beyond this; 0 is the full field (0)
%     'watertight_check'  'off', 'warn' or 'error' ('warn')
%   R is a struct: mesh (carrying sdf:distance), dims, origin and spacing
%   (1 x 3 rows), max_depth, num_banded and quality (see
%   mio.surface_watertight_check).

    o = mio.internal.opts(varargin, {'structure', 'voxel', 'resolution', [], 'cell_size', 0, ...
        'bounds', [], 'padding', 0, 'padding_relative', 0.1, 'root_resolution', 8, ...
        'max_depth', 4, 'band_cells', 1, 'record_levels', true, 'max_cells', 20000000, ...
        'sign', 'pseudonormal', 'location', 'corner', 'band', 0, ...
        'watertight_check', 'warn'}, 'mio.compute_sdf');
    structure = mio.internal.od_code(o.structure, {'voxel', 'octree'}, 'structure');
    sign = mio.internal.od_code(o.sign, {'unsigned', 'pseudonormal', 'winding-number'}, 'sign');
    location = mio.internal.od_sdf_location(o.location);
    check = mio.internal.od_code(o.watertight_check, {'off', 'warn', 'error'}, 'watertight check');
    r = meshioplusplus_mex('mio_compute_sdf', mesh, structure, double(o.resolution), ...
        double(o.cell_size), double(o.bounds), double(o.padding), ...
        double(o.padding_relative), double(o.root_resolution), double(o.max_depth), ...
        double(o.band_cells), logical(o.record_levels), double(o.max_cells), sign, ...
        location, double(o.band), check);
end
