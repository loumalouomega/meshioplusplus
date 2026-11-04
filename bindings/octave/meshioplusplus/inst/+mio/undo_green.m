function r = undo_green(coarse, fine)
%UNDO_GREEN  Undo the green closure of a refinement.
%   R = mio.undo_green(COARSE, FINE) restores, in FINE (a mio.refine result
%   carrying the 'refine:cell_id' hierarchy), every green-closure group back
%   to its COARSE parent verbatim. Fails by name when FINE has no hierarchy.
%   R is a struct: mesh, num_groups_undone and num_cells_removed.

    r = meshioplusplus_mex('mio_undo_green', coarse, fine);
end
