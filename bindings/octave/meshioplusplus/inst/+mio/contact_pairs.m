function r = contact_pairs(slave_mesh, slave_points, master_cells, varargin)
%CONTACT_PAIRS  Project a Point region onto a master Cell-region boundary.
%   R = mio.contact_pairs(SLAVE_MESH, SLAVE_POINTS, MASTER_CELLS, ...)
%   projects every point of the Point region SLAVE_POINTS to the closest
%   facet of the Cell region MASTER_CELLS. Options: 'master_mesh' ([]: the
%   slave mesh), 'tolerance' (0 derives a scale) and 'require_complete'
%   (false). Returns a struct of copies: slave_point and master_cell
%   (1-based columns), master_facet and master_subfacet (0-based columns),
%   local_coordinates, closest_point and normal (num_pairs x 3), gap
%   (column) and unmatched (1-based column). See doc/region_adjacency.md.

    o = mio.internal.opts(varargin, {'master_mesh', [], 'tolerance', 0, ...
        'require_complete', false}, 'mio.contact_pairs');
    r = meshioplusplus_mex('mio_contact_pairs', slave_mesh, ...
        mio.internal.str(slave_points, 'slave_points'), o.master_mesh, ...
        mio.internal.str(master_cells, 'master_cells'), double(o.tolerance), ...
        logical(o.require_complete));
end
