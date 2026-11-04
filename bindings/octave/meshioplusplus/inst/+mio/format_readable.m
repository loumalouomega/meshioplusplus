function tf = format_readable(format)
%FORMAT_READABLE  True when FORMAT can be read.

    tf = meshioplusplus_mex('mio_format_readable', mio.internal.str(format, 'format'));
end
