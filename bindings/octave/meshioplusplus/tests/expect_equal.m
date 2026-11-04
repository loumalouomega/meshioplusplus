function expect_equal(actual, expected, tol)
%EXPECT_EQUAL  Assert that ACTUAL equals EXPECTED (numeric values compared
%   as double within TOL, default 1e-12 relative; everything else with
%   isequal). Sizes must match exactly.
    if nargin < 3
        tol = 1e-12;
    end
    if isnumeric(actual) && isnumeric(expected) || islogical(actual) && islogical(expected)
        if ~isequal(size(actual), size(expected))
            error('size mismatch: got %s, expected %s', mat2str(size(actual)), ...
                  mat2str(size(expected)));
        end
        a = double(actual(:));
        e = double(expected(:));
        same = (a == e) | (isnan(a) & isnan(e)) | abs(a - e) <= tol * max(1, abs(e));
        if ~all(same)
            k = find(~same, 1);
            error('value mismatch at element %d: got %.17g, expected %.17g', k, a(k), e(k));
        end
    elseif ~isequal(actual, expected)
        error('values differ:\n  got:      %s\n  expected: %s', disp_str(actual), disp_str(expected));
    end
end

function s = disp_str(x)
    if ischar(x)
        s = ['''' x ''''];
    else
        s = strtrim(evalc('disp(x)'));
    end
end
