classdef MeshHandle < mio.Handle
%MESHHANDLE  A meshio++ mesh, owned by the native library.
%   Created by mio.mesh, mio.read and every operation that returns a mesh.
%   Use the mio.* functions on it (mio.points, mio.connectivity,
%   mio.point_data, ...), or copy it whole into a plain struct with TOSTRUCT
%   and back with mio.MeshHandle.fromStruct.

    methods
        function obj = MeshHandle(id)
            obj@mio.Handle(id);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  %s array of mio.MeshHandle\n', mat2str(size(obj)));
                return;
            end
            if ~obj.is_open()
                fprintf('  <mio.MeshHandle: released>\n');
            else
                fprintf('  <mio.MeshHandle: %g points, %g cells in %g block(s)>\n', ...
                        mio.num_points(obj), mio.num_cells(obj), mio.num_cell_blocks(obj));
            end
        end

        function s = toStruct(obj)
            %TOSTRUCT  Copy the mesh into a plain struct.
            %   S.points       num_points x dim double
            %   S.cells        struct array: type, data (num_cells x
            %                  nodes_per_cell, 1-based) or, for a ragged
            %                  block, a cell array of rows (polygons) or of
            %                  cells of face rows (polyhedra)
            %   S.point_data   struct array: name, data
            %   S.cell_data    struct array: name, data (cell array, one
            %                  array per cell block)
            %   S.field_data   struct array: name, data
            %   S.regions      struct array as returned by mio.regions
            s = mio.mesh_to_struct(obj);
        end
    end

    methods (Static)
        function m = fromStruct(s)
            %FROMSTRUCT  Build a new mesh from a struct shaped like TOSTRUCT's.
            m = mio.mesh_from_struct(s);
        end
    end
end
