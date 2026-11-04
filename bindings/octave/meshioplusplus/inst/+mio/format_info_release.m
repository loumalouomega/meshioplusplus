function format_info_release(info)
%FORMAT_INFO_RELEASE  Free a format side channel now.
%   mio.format_info_release(INFO) frees the native object behind the
%   mio.FormatInfo INFO immediately (otherwise it goes when the last
%   reference is cleared). Idempotent; anything but a mio.FormatInfo is an
%   error.

    meshioplusplus_mex('mio_format_info_release', info);
end
