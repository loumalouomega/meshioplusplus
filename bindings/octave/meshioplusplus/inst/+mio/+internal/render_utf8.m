function s = render_utf8(bytes)
%RENDER_UTF8  Decode a uint8 row of UTF-8 bytes from the MEX into a char row.
%   Octave's char arrays are UTF-8 bytes already; MATLAB's are UTF-16, so
%   there the bytes are decoded.

    bytes = reshape(bytes, 1, []);
    if exist('OCTAVE_VERSION', 'builtin') ~= 0
        s = char(bytes);
    else
        s = native2unicode(bytes, 'UTF-8');
    end
end
