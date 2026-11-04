function mio_build(varargin)
%MIO_BUILD  Compile the meshio++ MEX gateway (meshioplusplus_mex).
%   mio_build() compiles src/*.c of this package against the installed
%   meshio++ C library and writes meshioplusplus_mex.<mexext> next to this
%   file. It works in GNU Octave (mkoctfile --mex) and in MATLAB (mex).
%
%   The C library is located like the R binding's configure script does:
%     1. MESHIOPLUSPLUS_CFLAGS / MESHIOPLUSPLUS_LIBS  (fully manual)
%     2. MESHIOPLUSPLUS_HOME                          (an install prefix)
%     3. pkg-config meshioplusplus
%   or pass the prefix directly: mio_build('prefix', '/opt/meshioplusplus').
%
%   Octave users can instead `pkg install` the package tarball, which runs
%   the same build through src/Makefile.

    prefix = '';
    for k = 1:2:numel(varargin)
        if strcmp(varargin{k}, 'prefix')
            prefix = varargin{k + 1};
        else
            error('meshioplusplus:argument', 'mio_build: unknown option ''%s''', varargin{k});
        end
    end

    here = fileparts(mfilename('fullpath'));
    src = fullfile(here, '..', 'src');
    if ~exist(fullfile(src, 'mio_mex.c'), 'file')
        src = fullfile(here, 'src'); % installed layout
    end
    files = dir(fullfile(src, '*.c'));
    sources = cellfun(@(f) fullfile(src, f), {files.name}, 'UniformOutput', false);

    [cflags, libs] = locate(prefix);
    out = fullfile(here, ['meshioplusplus_mex.' mexext()]);
    is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if is_octave
        args = [{'--mex', '-o', out}, split_flags(cflags), sources, split_flags(libs)];
        [output, status] = mkoctfile(args{:});
        if status ~= 0
            error('meshioplusplus:build', 'mkoctfile failed:\n%s', output);
        end
    else
        args = [{'-R2018a', '-output', out}, split_flags(cflags), sources, split_flags(libs)];
        mex(args{:});
    end
    fprintf('meshio++: built %s\n', out);
end

function [cflags, libs] = locate(prefix)
    cflags = getenv('MESHIOPLUSPLUS_CFLAGS');
    libs = getenv('MESHIOPLUSPLUS_LIBS');
    if ~isempty(cflags) || ~isempty(libs)
        return;
    end
    if isempty(prefix)
        prefix = getenv('MESHIOPLUSPLUS_HOME');
    end
    if ~isempty(prefix)
        cflags = ['-I' fullfile(prefix, 'include')];
        libs = ['-L' fullfile(prefix, 'lib') ' -lmeshioplusplus'];
        return;
    end
    [s1, c] = system('pkg-config --cflags meshioplusplus');
    [s2, l] = system('pkg-config --libs meshioplusplus');
    if s1 == 0 && s2 == 0
        cflags = strtrim(c);
        libs = strtrim(l);
        return;
    end
    error('meshioplusplus:build', ['meshio++ C library not found. Build and install it ' ...
        '(build/configure.sh --c-api --build; cmake --install <tree> --prefix <p>), then ' ...
        'set MESHIOPLUSPLUS_HOME=<p>, put <p>/lib/pkgconfig on PKG_CONFIG_PATH, or call ' ...
        'mio_build(''prefix'', <p>).']);
end

function parts = split_flags(flags)
    parts = strsplit(strtrim(flags));
    parts = parts(~cellfun(@isempty, parts));
end
