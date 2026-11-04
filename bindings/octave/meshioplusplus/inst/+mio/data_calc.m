function out = data_calc(mesh, expression, output_name, varargin)
%DATA_CALC  Evaluate an elementwise expression into a new data array.
%   OUT = mio.data_calc(MESH, EXPRESSION, OUTPUT_NAME) evaluates EXPRESSION
%   over the arrays at the location and stores it as OUTPUT_NAME on a new
%   mesh. The grammar accepts + - * /, unary minus, parentheses, numeric
%   literals, array names and the functions abs, sqrt, min, max and norm.
%   Nothing else is evaluated: there is no arbitrary-code path. Options:
%     'location'   'point' (default), 'cell' or 'field'
%     'overwrite'  allow replacing an existing array (false)

    o = mio.internal.opts(varargin, {'location', 'point', 'overwrite', false}, 'mio.data_calc');
    out = meshioplusplus_mex('mio_data_calc', mesh, mio.internal.str(expression, 'expression'), ...
        mio.internal.choice(o.location, {'point', 'cell', 'field'}, 'location'), ...
        mio.internal.str(output_name, 'output_name'), logical(o.overwrite));
end
