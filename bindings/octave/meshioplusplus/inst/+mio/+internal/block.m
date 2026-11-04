function b = block(block)
%BLOCK  Validate a 1-based cell-block index and return it 0-based.

    if ~isnumeric(block) || ~isscalar(block) || isnan(block) || block < 1 || block ~= fix(block)
        error('meshioplusplus:argument', '`block` must be a single 1-based cell-block index');
    end
    b = double(block) - 1;
end
