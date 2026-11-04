function m = read(path, varargin)
%READ  Read a mesh file.
%   M = mio.read(PATH) infers the format from the extension (compressed
%   .gz/.bz2/.xz/.zst/.zip files are read transparently).
%   M = mio.read(PATH, FORMAT, ...) or name/value options:
%     'format'        explicit format name ([] to infer)
%     'points_only'   read the points only (false)
%     'metadata_only' read counts and names only (false)
%     'arrays'        cellstr of data arrays to read; [] reads every array,
%                     {} reads none
%     'mmap'          'auto', 'on' or 'off'
%     'time_step'     0-based time step of a transient file (0)
%     'lenient'       accept recoverable format violations (false)
%     'piece'         keep one 0-based piece of a partitioned file ([] merges)
%     'drop_ghosts'   drop ghost (halo) cells of .pvtu/.pvtp/.pvd (false)

    o = mio.internal.opts(varargin, {'format', [], 'points_only', false, ...
        'metadata_only', false, 'arrays', [], 'mmap', 'auto', 'time_step', 0, ...
        'lenient', false, 'piece', [], 'drop_ghosts', false}, 'mio.read');
    m = meshioplusplus_mex('mio_read', mio.internal.str(path, 'path'), ...
        mio.internal.str(o.format, 'format'), logical(o.points_only), ...
        logical(o.metadata_only), mio.internal.cellstr_arg(o.arrays, 'arrays'), ...
        mio.internal.choice(o.mmap, {'auto', 'on', 'off'}, 'mmap'), double(o.time_step), ...
        logical(o.lenient), o.piece, logical(o.drop_ghosts));
end
