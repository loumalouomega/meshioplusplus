function tf = format_info_is_open(info)
%FORMAT_INFO_IS_OPEN  True while a format side channel handle is live.
%   TF = mio.format_info_is_open(INFO) is false once INFO was released, and
%   for anything that is not a mio.FormatInfo.

    tf = meshioplusplus_mex('mio_format_info_is_open', info);
end
