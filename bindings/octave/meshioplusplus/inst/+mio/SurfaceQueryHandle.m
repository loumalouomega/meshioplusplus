classdef SurfaceQueryHandle < mio.Handle
%SURFACEQUERYHANDLE  A reusable surface query (mio.surface_query).

    methods
        function obj = SurfaceQueryHandle(id)
            obj@mio.Handle(id);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  %s array of mio.SurfaceQueryHandle\n', mat2str(size(obj)));
            elseif obj.is_open()
                fprintf('  <mio.SurfaceQueryHandle>\n');
            else
                fprintf('  <mio.SurfaceQueryHandle: released>\n');
            end
        end
    end
end
