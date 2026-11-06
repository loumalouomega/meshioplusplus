function r = remesh_volume(mesh, varargin)
%REMESH_VOLUME  Retetrahedralize a volume (or closed surface) by isosurface stuffing.
%   R = mio.remesh_volume(MESH, ...) -- mio.remesh's volumetric sibling.
%   Exactly one of 'resolution' / 'cell_size' must be given, sizing a
%   body-centered cubic (BCC) lattice whose uncut tets have dihedral angles
%   from a fixed, mesh-size-independent set. 'warp_fraction' moves lattice
%   vertices near the surface onto it, trading a small, measured chance of
%   non-manifold boundary edges (num_non_manifold_edges) for substantially
%   better boundary tet quality; 0 disables warping and gives an exactly
%   watertight but lower-quality boundary. See doc/remesh_volume.md.
%
%   Attribution: implemented from the PUBLISHED DESCRIPTION of Labelle &
%   Shewchuk, "Isosurface Stuffing" (SIGGRAPH 2007) only -- neither the
%   paper's reference implementation (Stellar) nor TetGen (AGPL-3.0) is read
%   or vendored here.
%
%   Options (defaults in parentheses):
%     'resolution'        root lattice cell counts [nx ny nz] ([])
%     'cell_size'         cubic cell size of the root lattice (0, unset)
%     'bounds'            explicit [xlo ylo zlo xhi yhi zhi]; [] uses the
%                         surface's own bounding box
%     'padding'           padding on every side, world units (0)
%     'padding_relative'  padding as a fraction of the diagonal (0.1)
%     'max_cells'         refuse a root lattice above this many cells (20000000)
%     'max_tets'          refuse an output with more tets than this, checked
%                         after cutting (20000000)
%     'warp_fraction'     fraction of a vertex's shortest incident edge within
%                         which it may be warped onto the surface (0.35)
%     'sign'              'pseudonormal' or 'winding-number' ('pseudonormal')
%     'watertight_check'  'off', 'warn' or 'error' ('warn')
%     'surface_region'    restrict the surface to a named Cell region ('' all)
%     'grid_cell_size'    ignored (sized the search grid a BVH replaced) (0)
%     'max_winding_work'  refuse a 'winding-number' query above this cost (2e9)
%   R is a struct: mesh, num_tets, num_vertices_warped, num_tets_rejected,
%   num_non_manifold_edges and input_quality (the verdict for the INPUT
%   surface, not the output).

    o = mio.internal.opts(varargin, {'resolution', [], 'cell_size', 0, 'bounds', [], ...
        'padding', 0, 'padding_relative', 0.1, 'max_cells', 20000000, ...
        'max_tets', 20000000, 'warp_fraction', 0.35, 'sign', 'pseudonormal', ...
        'watertight_check', 'warn', 'surface_region', '', 'grid_cell_size', 0, ...
        'max_winding_work', 2e9}, 'mio.remesh_volume');
    sign = mio.internal.od_code(o.sign, {'unsigned', 'pseudonormal', 'winding-number'}, 'sign');
    check = mio.internal.od_code(o.watertight_check, {'off', 'warn', 'error'}, 'watertight check');
    r = meshioplusplus_mex('mio_remesh_volume', mesh, double(o.resolution), ...
        double(o.cell_size), double(o.bounds), double(o.padding), ...
        double(o.padding_relative), double(o.max_cells), double(o.max_tets), ...
        double(o.warp_fraction), sign, check, mio.internal.str(o.surface_region, ...
        'surface_region'), double(o.grid_cell_size), double(o.max_winding_work));
end
