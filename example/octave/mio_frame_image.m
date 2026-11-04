function f = mio_frame_image(mesh, varargin)
%MIO_FRAME_IMAGE  Draw MESH with the software rasterizer and show it inline.
%   F = mio_frame_image(MESH, 'name', value, ...) passes the options of
%   mio.render ('view', 'color_by', 'streamlines', 'cutaway', ...) through,
%   with 'width' 560 and 'height' 420 unless given, shows the frame with
%   image() in a new figure the size of the frame, and returns it. The
%   kernel captures the figure as an inline PNG. Notebook helper, not part
%   of the package API.

    width = 560;
    height = 420;
    rest = {};
    for k = 1:2:numel(varargin)
        switch varargin{k}
            case 'width'
                width = varargin{k + 1};
            case 'height'
                height = varargin{k + 1};
            otherwise
                rest(end + 1:end + 2) = varargin(k:k + 1); %#ok<AGROW>
        end
    end
    f = mio.render(mesh, 'width', width, 'height', height, rest{:});
    figure('position', [100 100 width height]);
    axes('position', [0 0 1 1]);
    image(f.rgba(:, :, 1:3));
    axis image off
end
