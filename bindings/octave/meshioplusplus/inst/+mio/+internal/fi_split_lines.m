function c = fi_split_lines(s)
%FI_SPLIT_LINES  A newline-separated string as an n x 1 cellstr ('' -> {}).

    if isempty(s)
        c = cell(0, 1);
        return;
    end
    c = strsplit(s, sprintf('\n'));
    c = c(:);
end
