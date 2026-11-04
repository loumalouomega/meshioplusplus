function tf = has_cutcells()
%HAS_CUTCELLS  Whether the library carries CutCells.
%   mio.cut_level_set and mio.cut_quadrature need it (a build with
%   MESHIOPLUSPLUS_WITH_CUTCELLS).

    tf = meshioplusplus_mex('mio_has_cutcells');
end
