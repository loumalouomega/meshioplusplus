function o = opts(args, spec, fname)
%OPTS  Parse the optional arguments of a mio.* function.
%   O = mio.internal.opts(ARGS, SPEC, FNAME) parses the cell ARGS (a
%   function's VARARGIN) against SPEC, a cell {'name1', default1, 'name2',
%   default2, ...} listing the optional arguments in the order the R binding
%   declares them, and returns a struct with one field per name.
%
%   Optional arguments may be given positionally, in SPEC order, and then as
%   name/value pairs -- the R binding's calling convention:
%
%       mio.connectivity(m, 2)
%       mio.read('a.msh', 'gmsh')
%       mio.read('a.msh', 'format', 'gmsh', 'lenient', true)
%
%   A char argument that names an option starts the name/value part, so a
%   positional string that happens to equal an option name must be passed as
%   a pair instead. Unknown names, a missing value and too many positional
%   arguments are errors (identifier meshioplusplus:argument).

    if nargin < 3
        fname = 'mio';
    end
    names = spec(1:2:end);
    o = struct();
    for k = 1:numel(names)
        o.(names{k}) = spec{2 * k};
    end

    n = numel(args);
    i = 1;
    pos = 0;
    while i <= n
        a = args{i};
        if ischar(a) && isrow(a) && any(strcmp(a, names)) && i < n
            break;
        end
        if isstruct(a) && isscalar(a) && i == n && n > 0 && pos == 0 && ~any(strcmp('struct', names))
            % A single options struct: every field must be a known option.
            f = fieldnames(a);
            for k = 1:numel(f)
                if ~any(strcmp(f{k}, names))
                    error('meshioplusplus:argument', '%s: unknown option ''%s''', fname, f{k});
                end
                o.(f{k}) = a.(f{k});
            end
            return;
        end
        pos = pos + 1;
        if pos > numel(names)
            % The likeliest cause is a misspelt option name taken as a value.
            for k = i - 1:-1:1
                if ischar(args{k})
                    error('meshioplusplus:argument', ...
                          '%s: too many arguments (unknown option ''%s''?)', fname, args{k});
                end
            end
            error('meshioplusplus:argument', '%s: too many arguments', fname);
        end
        o.(names{pos}) = a;
        i = i + 1;
    end
    if mod(n - i + 1, 2) ~= 0
        error('meshioplusplus:argument', '%s: options must be name/value pairs', fname);
    end
    while i <= n
        name = args{i};
        if ~ischar(name) || ~any(strcmp(name, names))
            if ischar(name)
                error('meshioplusplus:argument', '%s: unknown option ''%s''', fname, name);
            end
            error('meshioplusplus:argument', '%s: option names must be strings', fname);
        end
        o.(name) = args{i + 1};
        i = i + 2;
    end
end
