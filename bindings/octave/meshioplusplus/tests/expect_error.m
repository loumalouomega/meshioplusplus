function expect_error(fn, pattern)
%EXPECT_ERROR  Assert that FN() raises, optionally with a message matching
%   the regular expression PATTERN.
    try
        fn();
    catch err
        if nargin > 1 && isempty(regexp(err.message, pattern, 'once'))
            error('expected an error matching ''%s'', got: %s', pattern, err.message);
        end
        return;
    end
    error('expected an error, but the call succeeded');
end
