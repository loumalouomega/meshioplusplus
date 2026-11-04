classdef XdmfSeries < mio.Handle
%XDMFSERIES  A transient XDMF writer (mio.xdmf_series).

    methods
        function obj = XdmfSeries(id)
            obj@mio.Handle(id);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  %s array of mio.XdmfSeries\n', mat2str(size(obj)));
            elseif obj.is_open()
                fprintf('  <mio.XdmfSeries>\n');
            else
                fprintf('  <mio.XdmfSeries: released>\n');
            end
        end
    end
end
