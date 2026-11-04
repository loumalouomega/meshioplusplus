function tf = reader_supports_options(format)
%READER_SUPPORTS_OPTIONS  True when the reader of FORMAT honours selective-read options.

    tf = meshioplusplus_mex('mio_reader_supports_options', mio.internal.str(format, 'format'));
end
