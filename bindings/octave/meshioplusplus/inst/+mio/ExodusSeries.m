classdef ExodusSeries < mio.Handle
%EXODUSSERIES  A stateful Exodus II series writer (mio.exodus_series).
%   Created by mio.exodus_series; see that function for the lifecycle. Displays as
%   <mio.ExodusSeries: N step(s) written> (", finalized" once finalized), or
%   <mio.ExodusSeries: released>.

    methods
        function obj = ExodusSeries(id)
            obj@mio.Handle(id);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  %s array of mio.ExodusSeries\n', mat2str(size(obj)));
            elseif ~obj.is_open()
                fprintf('  <mio.ExodusSeries: released>\n');
            else
                n = mio.exodus_series_num_steps(obj);
                if mio.exodus_series_finalized(obj)
                    state = ', finalized';
                else
                    state = '';
                end
                fprintf('  <mio.ExodusSeries: %g step(s) written%s>\n', n, state);
            end
        end
    end
end
