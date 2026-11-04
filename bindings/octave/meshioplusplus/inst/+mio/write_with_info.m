function write_with_info(mesh, info, path, varargin)
%WRITE_WITH_INFO  Write a mesh together with a format side channel.
%   mio.write_with_info(MESH, INFO, PATH) writes MESH to PATH and puts back
%   the side channel INFO (a mio.FormatInfo from mio.read_with_info,
%   mio.iga_refine or mio.iga_sample). With INFO = [] it is exactly
%   mio.write. Option: 'format' (explicit format name, [] to infer).

    o = mio.internal.opts(varargin, {'format', []}, 'mio.write_with_info');
    meshioplusplus_mex('mio_write_with_info', mesh, info, mio.internal.str(path, 'path'), ...
        mio.internal.str(o.format, 'format'));
end
