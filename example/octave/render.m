function f = render(mesh, varargin)
%RENDER  Show MESH inline: smooth shading, feature edges, white background.
%   render(MESH) is the counterpart of the R notebooks' render(): where R
%   goes through the SVG writer, whose flat-API entry point takes no styling,
%   this goes through mio.render, which takes every option. Extra
%   name/value pairs ('color_by', 'cmap', 'colorbar', 'view', ...) are passed
%   to mio.render. Notebook helper, not part of the package API.

    f = mio_frame_image(mesh, 'width', 480, 'height', 360, 'shading', 'smooth', ...
        'edges', 'feature', 'background', [255 255 255], varargin{:});
    if nargout == 0
        clear f
    end
end
