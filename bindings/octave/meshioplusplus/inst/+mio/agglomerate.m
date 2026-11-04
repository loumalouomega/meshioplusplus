function r = agglomerate(mesh, varargin)
%AGGLOMERATE  Coarsen a mesh by merging cells into polyhedral groups.
%   R = mio.agglomerate(MESH, ...). Options:
%     'target_group_size'    approximate member cells per output group, at
%                            least 1 (1 groups every cell by itself) (8)
%     'merge_coplanar_faces' fuse coplanar shared and boundary faces into
%                            one polygon (false)
%     'coplanar_angle'       largest normal deviation (degrees) still coplanar (1)
%     'min_sphericity'       refuse an absorption that would drop the
%                            group's sphericity below this; 0 disables (0)
%   R is a struct: mesh, cell_map (one 1-based entry per input cell: its
%   output cell), num_faces_merged and num_rejected.

    o = mio.internal.opts(varargin, {'target_group_size', 8, 'merge_coplanar_faces', false, ...
        'coplanar_angle', 1, 'min_sphericity', 0}, 'mio.agglomerate');
    r = meshioplusplus_mex('mio_agglomerate', mesh, double(o.target_group_size), ...
        logical(o.merge_coplanar_faces), double(o.coplanar_angle), double(o.min_sphericity));
end
