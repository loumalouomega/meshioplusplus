function k = choice(value, choices, what)
%CHOICE  0-based index of VALUE in the cellstr CHOICES (R's match.arg).
%   An exact match is required; the error lists the accepted values.

    if ischar(value)
        hit = find(strcmp(value, choices), 1);
        if ~isempty(hit)
            k = hit - 1;
            return;
        end
    end
    error('meshioplusplus:argument', '`%s` must be one of: %s', what, ...
          strjoin(choices, ', '));
end
