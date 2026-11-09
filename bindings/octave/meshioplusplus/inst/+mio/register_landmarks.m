function r = register_landmarks(source, target, varargin)
%REGISTER_LANDMARKS  The rigid motion laying picked points onto their targets.
%   R = mio.register_landmarks(SOURCE, TARGET) solves the weighted Horn
%   problem for n x 3 matrices of picked pairs (n >= 3, the sources not
%   collinear; roadmap 7.3.2). Option 'weights' (n non-negative weights, []).
%   Returns a struct: matrix (4x4), residuals (n x 1) and rms.

    o = mio.internal.opts(varargin, {'weights', []}, 'mio.register_landmarks');
    w = [];
    if ~isempty(o.weights)
        w = double(o.weights(:));
    end
    r = meshioplusplus_mex('mio_register_landmarks', double(source), double(target), w);
end
