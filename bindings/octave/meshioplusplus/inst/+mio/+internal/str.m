function s = str(value, what)
%STR  A string argument as a char row vector ([] stays [] -- the R NULL).
%   Accepts char and (in MATLAB) string scalars.

    if isempty(value) && ~ischar(value)
        s = [];
        return;
    end
    if isa(value, 'string')
        value = char(value);
    end
    if ~ischar(value)
        if nargin < 2
            what = 'argument';
        end
        error('meshioplusplus:argument', '`%s` must be a string', what);
    end
    s = value;
end
