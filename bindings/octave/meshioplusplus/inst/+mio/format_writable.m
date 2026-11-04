function tf = format_writable(format)
%FORMAT_WRITABLE  True when FORMAT can be written.

    tf = meshioplusplus_mex('mio_format_writable', mio.internal.str(format, 'format'));
end
