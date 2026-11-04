function m = mesh_from_struct(s)
%MESH_FROM_STRUCT  Build a new mesh from a struct shaped like
%   mio.MeshHandle/toStruct's output. Only 'points' is required; 'cells',
%   'point_data', 'cell_data', 'field_data' and 'regions' are optional. A
%   cell block's data is a connectivity matrix, a cell array of rows
%   (polygons), or a cell array of cells of face rows (polyhedra).

    if ~isstruct(s) || ~isscalar(s) || ~isfield(s, 'points')
        error('meshioplusplus:argument', 'expected a struct with at least a ''points'' field');
    end
    m = mio.mesh();
    mio.set_points(m, s.points);
    if isfield(s, 'cells')
        for b = 1:numel(s.cells)
            data = s.cells(b).data;
            if iscell(data)
                if ~isempty(data) && iscell(data{1})
                    mio.add_polyhedron_block(m, s.cells(b).type, data);
                else
                    mio.add_polygon_block(m, s.cells(b).type, data);
                end
            else
                mio.add_cell_block(m, s.cells(b).type, data);
            end
        end
    end
    if isfield(s, 'point_data')
        for k = 1:numel(s.point_data)
            mio.add_point_data(m, s.point_data(k).name, s.point_data(k).data);
        end
    end
    if isfield(s, 'cell_data')
        for k = 1:numel(s.cell_data)
            data = s.cell_data(k).data;
            if ~iscell(data)
                data = {data};
            end
            for b = 1:numel(data)
                mio.append_cell_data(m, s.cell_data(k).name, data{b});
            end
        end
    end
    if isfield(s, 'field_data')
        for k = 1:numel(s.field_data)
            mio.add_field_data(m, s.field_data(k).name, s.field_data(k).data);
        end
    end
    if isfield(s, 'regions')
        for k = 1:numel(s.regions)
            r = s.regions(k);
            mio.add_region(m, r.name, r.kind, r.entries, 'dim', r.dim, 'tag', r.tag);
        end
    end
end
