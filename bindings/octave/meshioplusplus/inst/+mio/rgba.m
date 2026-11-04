function c = rgba(r, g, b, varargin)
%RGBA  A colour as the 0xRRGGBBAA number the render options take.
%   C = mio.rgba(R, G, B) or mio.rgba(R, G, B, A) with channels in 0..255;
%   A defaults to 255 (opaque). C is a double: it can exceed the int32
%   range. Inputs may be arrays of equal size.

    o = mio.internal.opts(varargin, {'a', 255}, 'mio.rgba');
    c = ((double(r) * 256 + double(g)) * 256 + double(b)) * 256 + double(o.a);
end
