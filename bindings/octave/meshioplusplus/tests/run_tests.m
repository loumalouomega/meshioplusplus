function failures = run_tests(pattern)
%RUN_TESTS  Run the meshio++ Octave/MATLAB test suite.
%   RUN_TESTS() runs every tests/test_*.m; RUN_TESTS('render') only the
%   files whose name contains 'render'. Each test file is a function
%   returning a cell array of {name, @fn} pairs (see test_core.m). A test
%   passes when its function returns, is skipped when it raises
%   meshioplusplus:skip (via skip_if_not) and fails on any other error.
%   Returns the number of failures, so CI runs
%       octave --no-gui --eval "exit(run_tests() > 0)"
%       matlab -batch "exit(run_tests() > 0)"

    if nargin < 1
        pattern = '';
    end
    here = fileparts(mfilename('fullpath'));
    addpath(here);
    inst = fullfile(here, '..', 'inst');
    if exist(fullfile(inst, '+mio'), 'dir') && isempty(which('mio.read'))
        addpath(inst);
    end
    files = dir(fullfile(here, 'test_*.m'));
    passed = 0;
    skipped = 0;
    failures = 0;
    for f = 1:numel(files)
        [~, name] = fileparts(files(f).name);
        if ~isempty(pattern) && isempty(strfind(name, pattern))
            continue;
        end
        tests = feval(name);
        for t = 1:size(tests, 1)
            label = sprintf('%s: %s', name, tests{t, 1});
            try
                tests{t, 2}();
                passed = passed + 1;
                fprintf('  ok    %s\n', label);
            catch err
                if strcmp(err.identifier, 'meshioplusplus:skip')
                    skipped = skipped + 1;
                    fprintf('  skip  %s (%s)\n', label, err.message);
                else
                    failures = failures + 1;
                    fprintf('  FAIL  %s\n        %s\n', label, err.message);
                    for s = 1:min(numel(err.stack), 3)
                        fprintf('        at %s:%d\n', err.stack(s).name, err.stack(s).line);
                    end
                end
            end
        end
    end
    fprintf('\n%d passed, %d skipped, %d failed\n', passed, skipped, failures);
end
