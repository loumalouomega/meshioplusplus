function s = mesh_to_struct(mesh)
%MESH_TO_STRUCT  Copy MESH into a plain struct (see mio.MeshHandle/toStruct).
%   Data names such as 'gmsh:physical' are not valid struct field names, so
%   the data sections are struct arrays of name/data pairs rather than
%   structs keyed by name.

    s = struct();
    s.points = mio.points(mesh);
    nb = mio.num_cell_blocks(mesh);
    cells = struct('type', cell(nb, 1), 'data', cell(nb, 1));
    for b = 1:nb
        info = mio.cell_block_info(mesh, b);
        cells(b).type = mio.cell_block_type(mesh, b);
        if info.is_polyhedron
            cells(b).data = mio.polyhedron_block(mesh, b);
        elseif info.is_ragged
            cells(b).data = mio.polygon_block(mesh, b);
        else
            cells(b).data = mio.connectivity(mesh, b);
        end
    end
    s.cells = cells;

    names = mio.point_data_names(mesh);
    s.point_data = struct('name', names, 'data', cell(size(names)));
    for k = 1:numel(names)
        s.point_data(k).data = mio.point_data(mesh, names{k});
    end

    names = mio.cell_data_names(mesh);
    s.cell_data = struct('name', names, 'data', cell(size(names)));
    for k = 1:numel(names)
        nblk = mio.cell_data_num_blocks(mesh, names{k});
        per_block = cell(nblk, 1);
        for b = 1:nblk
            per_block{b} = mio.cell_data(mesh, names{k}, b);
        end
        s.cell_data(k).data = per_block;
    end

    names = mio.field_data_names(mesh);
    s.field_data = struct('name', names, 'data', cell(size(names)));
    for k = 1:numel(names)
        s.field_data(k).data = mio.field_data(mesh, names{k});
    end

    s.regions = mio.regions(mesh);
end
