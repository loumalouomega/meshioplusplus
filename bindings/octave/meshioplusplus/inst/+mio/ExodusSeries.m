classdef ExodusSeries < mio.Handle
%EXODUSSERIES  A stateful Exodus II series writer (mio.exodus_series).

    methods
        function obj = ExodusSeries(id)
            obj@mio.Handle(id);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  %s array of mio.ExodusSeries\n', mat2str(size(obj)));
            elseif obj.is_open()
                fprintf('  <mio.ExodusSeries>\n');
            else
                fprintf('  <mio.ExodusSeries: released>\n');
            end
        end
    end
end
