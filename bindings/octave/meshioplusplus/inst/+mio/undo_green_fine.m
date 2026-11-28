function r = undo_green_fine(fine)
%UNDO_GREEN_FINE  Undo the green closures of a refined mesh from its own hierarchy.
%   R = mio.undo_green_fine(FINE) is the single-mesh form of mio.undo_green
%   (roadmap 9.6.1): the green groups come from the refine:cell_id and
%   refine:parent_id arrays FINE carries (a 'record_hierarchy', true refine).
%   Returns a struct: mesh, num_groups_undone, num_cells_removed.

    r = meshioplusplus_mex('mio_undo_green_fine', fine);
end
