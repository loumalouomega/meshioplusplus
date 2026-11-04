classdef FormatInfo < mio.Handle
%FORMATINFO  A format side channel (Gmsh periodic info, NURBS / B-rep / MDPA data).

    methods
        function obj = FormatInfo(id)
            obj@mio.Handle(id);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  %s array of mio.FormatInfo\n', mat2str(size(obj)));
            elseif obj.is_open()
                fprintf('  <mio.FormatInfo>\n');
            else
                fprintf('  <mio.FormatInfo: released>\n');
            end
        end
    end
end
