function f = format_info_format(info)
%FORMAT_INFO_FORMAT  The format name a side channel belongs to.
%   F = mio.format_info_format(INFO) returns e.g. 'mdpa', 'gmsh', 'mfem' or
%   'cad_json' for a mio.FormatInfo.

    f = meshioplusplus_mex('mio_format_info_format', info);
end
