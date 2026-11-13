function r = detect(path)
%DETECT  A ranked format detection report for PATH.
%   R = mio.detect(PATH) is a struct array with fields format, confidence
%   ('certain': a content signature matched; 'likely': the extension names it and
%   its metadata reads; 'guess': the extension alone) and reason, best first;
%   empty when nothing matches.

    text = meshioplusplus_mex('mio_detect', mio.internal.str(path, 'path'));
    r = struct('format', {}, 'confidence', {}, 'reason', {});
    lines = strsplit(text, sprintf('\n'));
    for k = 1:numel(lines)
        if isempty(lines{k})
            continue;
        end
        f = strsplit(lines{k}, sprintf('\t'));
        if numel(f) < 3
            continue;
        end
        r(end + 1) = struct('format', f{1}, 'confidence', f{2}, ...
                            'reason', strjoin(f(3:end), sprintf('\t')));
    end
end
