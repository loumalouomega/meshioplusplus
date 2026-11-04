function out = data_condition(mesh, location, varargin)
%DATA_CONDITION  Clamp, normalize or standardize data arrays.
%   OUT = mio.data_condition(MESH, LOCATION) conditions the arrays at
%   LOCATION ('point', 'cell' or 'field') on a new mesh. Options:
%     'names'            arrays to condition ([] or {} for every array)
%     'cond_mode'        'clamp' (default), 'normalize' or 'standardize'
%     'lo', 'hi'         clamp bounds, or the normalize target range (0, 1)
%     'scope'            'component' (each trailing component independently,
%                        default) or 'magnitude' (rescale whole rows,
%                        preserving direction)
%     'nan_policy'       'ignore' (default), 'replace' or 'fail'; non-finite
%                        values are always excluded from reductions
%     'nan_replacement'  used when nan_policy is 'replace' (0)
%     'suffix'           write to name + suffix instead of in place ([])

    o = mio.internal.opts(varargin, {'names', [], 'cond_mode', 'clamp', 'lo', 0, 'hi', 1, ...
        'scope', 'component', 'nan_policy', 'ignore', 'nan_replacement', 0, 'suffix', []}, ...
        'mio.data_condition');
    md = mio.internal.choice(o.cond_mode, {'clamp', 'normalize', 'standardize'}, 'cond_mode');
    sc = mio.internal.choice(o.scope, {'component', 'magnitude'}, 'scope');
    np = mio.internal.choice(o.nan_policy, {'ignore', 'replace', 'fail'}, 'nan_policy');
    out = meshioplusplus_mex('mio_data_condition', mesh, ...
        mio.internal.choice(location, {'point', 'cell', 'field'}, 'location'), ...
        mio.internal.cellstr_arg(o.names, 'names'), md, double(o.lo), double(o.hi), sc, np, ...
        double(o.nan_replacement), mio.internal.str(o.suffix, 'suffix'));
end
