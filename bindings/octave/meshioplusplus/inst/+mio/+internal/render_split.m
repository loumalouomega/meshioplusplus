function [o, pairs] = render_split(args, spec, fname, positional)
%RENDER_SPLIT  Split a render wrapper's VARARGIN into its own options and
%   the render options.
%   [O, PAIRS] = mio.internal.render_split(ARGS, SPEC, FNAME, POSITIONAL)
%   parses ARGS like mio.internal.opts against SPEC ({'name', default, ...}),
%   except that a name/value pair whose name is not in SPEC is a render
%   option (R's `...`): those are returned in order as the 1 x 2k cell PAIRS
%   the MEX takes, and the MEX rejects an unknown name. With POSITIONAL
%   true, leading arguments fill SPEC in order until the first option name;
%   with it false (an R signature whose `...` comes first) every option is a
%   name/value pair. A single struct argument supplies every option by field.

    names = spec(1:2:end);
    o = struct();
    for k = 1:numel(names)
        o.(names{k}) = spec{2 * k};
    end
    pairs = cell(1, 0);
    n = numel(args);
    if n == 1 && isstruct(args{1}) && isscalar(args{1})
        s = args{1};
        f = fieldnames(s);
        args = cell(1, 2 * numel(f));
        for k = 1:numel(f)
            args{2 * k - 1} = f{k};
            args{2 * k} = s.(f{k});
        end
        n = numel(args);
    end

    i = 1;
    if positional
        known = [names, mio.internal.render_names()];
        pos = 0;
        while i <= n
            a = args{i};
            if ischar(a) && isrow(a) && any(strcmp(a, known))
                break;
            end
            pos = pos + 1;
            if pos > numel(names)
                error('meshioplusplus:argument', '%s: too many positional arguments', fname);
            end
            if ischar(a) && ~ischar(spec{2 * pos})
                error('meshioplusplus:argument', '%s: unknown option ''%s''', fname, a);
            end
            o.(names{pos}) = a;
            i = i + 1;
        end
    end
    if mod(n - i + 1, 2) ~= 0
        error('meshioplusplus:argument', '%s: options must be name/value pairs', fname);
    end
    while i <= n
        name = args{i};
        if ~ischar(name) || ~isrow(name)
            error('meshioplusplus:argument', '%s: option names must be strings', fname);
        end
        if any(strcmp(name, names))
            o.(name) = args{i + 1};
        else
            pairs{end + 1} = name; %#ok<AGROW>
            pairs{end + 1} = args{i + 1}; %#ok<AGROW>
        end
        i = i + 2;
    end
end
