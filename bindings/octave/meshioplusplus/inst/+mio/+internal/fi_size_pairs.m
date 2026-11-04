function [idx, sz] = fi_size_pairs(x, what)
%FI_SIZE_PAIRS  (1-based index, size) pairs of mio.mesh_cad_surface.
%   X is [] (none) or an (n x 2) matrix whose first column holds 1-based
%   face or edge indices and whose second the sizes -- the Octave form of
%   R's named vector / two-column matrix. Returns two n x 1 columns.

    if isempty(x)
        idx = zeros(0, 1);
        sz = zeros(0, 1);
        return;
    end
    if ~isnumeric(x) || ndims(x) ~= 2 || size(x, 2) ~= 2
        error('meshioplusplus:argument', ...
              '`%s` must be an (n x 2) matrix of [1-based index, size] rows', what);
    end
    idx = double(x(:, 1));
    sz = double(x(:, 2));
end
