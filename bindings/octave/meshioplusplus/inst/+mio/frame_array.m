function a = frame_array(frame)
%FRAME_ARRAY  A rendered frame as a double array in [0, 1].
%   A = mio.frame_array(F) takes a frame from mio.render and returns its
%   pixels as a height x width x 4 double array (red, green, blue, alpha),
%   row 1 on top: image(A(:,:,1:3)) shows it, as R's rasterImage() takes
%   the same array there.

    if ~isstruct(frame) || ~isfield(frame, 'rgba')
        error('meshioplusplus:argument', '`frame` must be a frame from mio.render');
    end
    a = double(frame.rgba) / 255;
end
