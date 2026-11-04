classdef SequenceHandle < mio.Handle
%SEQUENCEHANDLE  A transient / multi-file dataset (mio.sequence).

    methods
        function obj = SequenceHandle(id)
            obj@mio.Handle(id);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  %s array of mio.SequenceHandle\n', mat2str(size(obj)));
            elseif obj.is_open()
                fprintf('  <mio.SequenceHandle>\n');
            else
                fprintf('  <mio.SequenceHandle: released>\n');
            end
        end
    end
end
