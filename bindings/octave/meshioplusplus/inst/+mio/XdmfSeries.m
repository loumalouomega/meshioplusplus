classdef XdmfSeries < mio.Handle
%XDMFSERIES  A transient XDMF writer (mio.xdmf_series).
%   Created by mio.xdmf_series; see that function for the lifecycle. Displays as
%   <mio.XdmfSeries: N step(s) written> (", finalized" once finalized), or
%   <mio.XdmfSeries: released>.

    methods
        function obj = XdmfSeries(id)
            obj@mio.Handle(id);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  %s array of mio.XdmfSeries\n', mat2str(size(obj)));
            elseif ~obj.is_open()
                fprintf('  <mio.XdmfSeries: released>\n');
            else
                n = mio.xdmf_series_num_steps(obj);
                if mio.xdmf_series_finalized(obj)
                    state = ', finalized';
                else
                    state = '';
                end
                fprintf('  <mio.XdmfSeries: %g step(s) written%s>\n', n, state);
            end
        end
    end
end
