function r = split_warnings(r)
%SPLIT_WARNINGS  R.warnings, newline-separated text, as a column cellstr.

    if isempty(r.warnings)
        r.warnings = cell(0, 1);
    else
        r.warnings = strsplit(r.warnings, char(10));
        r.warnings = r.warnings(:);
    end
end
