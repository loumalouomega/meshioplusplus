function mesh = add_polyhedron_block(mesh, cell_type, cells)
%ADD_POLYHEDRON_BLOCK  Append a polyhedron block: CELLS is a cell array (one per cell) of
%   cell arrays of 1-based face node vectors. Returns MESH.

    meshioplusplus_mex('mio_add_polyhedron_block', mesh, mio.internal.str(cell_type, 'cell_type'), cells);
end
