function convert(in_path, out_path, varargin)
%CONVERT  Convert a mesh file without materialising it here.
%   Options: 'in_format', 'out_format' ([] to infer).

    o = mio.internal.opts(varargin, {'in_format', [], 'out_format', []}, 'mio.convert');
    meshioplusplus_mex('mio_convert', mio.internal.str(in_path, 'in_path'), ...
        mio.internal.str(o.in_format, 'in_format'), mio.internal.str(out_path, 'out_path'), ...
        mio.internal.str(o.out_format, 'out_format'));
end
