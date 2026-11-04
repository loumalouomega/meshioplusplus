function skip_if_not(condition, reason)
%SKIP_IF_NOT  Skip the current test unless CONDITION holds.
    if ~condition
        if nargin < 2
            reason = 'condition not met';
        end
        error('meshioplusplus:skip', '%s', reason);
    end
end
