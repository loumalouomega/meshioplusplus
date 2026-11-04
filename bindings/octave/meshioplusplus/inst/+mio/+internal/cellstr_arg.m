function c = cellstr_arg(value, what)
%CELLSTR_ARG  A list-of-strings argument as a cellstr. [] stays [] (the R
%   NULL); a single char row becomes a one-element cell; {} stays {}.

    if isempty(value) && isnumeric(value)
        c = [];
        return;
    end
    if isa(value, 'string')
        value = cellstr(value);
    end
    if ischar(value)
        c = {value};
        return;
    end
    if ~iscellstr(value)
        if nargin < 2
            what = 'argument';
        end
        error('meshioplusplus:argument', '`%s` must be a string or a cell array of strings', what);
    end
    c = value;
end
