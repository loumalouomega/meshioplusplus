function s = sniff_format(path)
%SNIFF_FORMAT  The format detected from PATH's content, or ''.

    s = meshioplusplus_mex('mio_sniff_format', mio.internal.str(path, 'path'));
end
