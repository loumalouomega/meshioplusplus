function code = od_code(value, names, what)
%OD_CODE  0-based index of VALUE in NAMES, raising "unknown WHAT 'VALUE'"
%   otherwise -- the R binding's `if (!x %in% names(table)) stop(...)`.

    hit = [];
    if ischar(value)
        hit = find(strcmp(value, names), 1);
    end
    if isempty(hit)
        if ischar(value)
            shown = value;
        elseif isnumeric(value) || islogical(value)
            shown = mat2str(value);
        else
            shown = class(value);
        end
        error('meshioplusplus:argument', 'unknown %s ''%s''', what, shown);
    end
    code = hit - 1;
end
