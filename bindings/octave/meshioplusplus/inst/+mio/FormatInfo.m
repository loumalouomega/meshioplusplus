classdef FormatInfo < mio.Handle
%FORMATINFO  A format side channel (Gmsh periodic info, NURBS / B-rep / MDPA data).
%   Made by mio.read_with_info, mio.iga_refine and mio.iga_sample. Displays
%   as <mio.FormatInfo: FORMAT>, or <mio.FormatInfo: released>.

    methods
        function obj = FormatInfo(id)
            obj@mio.Handle(id);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  %s array of mio.FormatInfo\n', mat2str(size(obj)));
            elseif obj.is_open()
                fprintf('  <mio.FormatInfo: %s>\n', mio.format_info_format(obj));
            else
                fprintf('  <mio.FormatInfo: released>\n');
            end
        end
    end
end
