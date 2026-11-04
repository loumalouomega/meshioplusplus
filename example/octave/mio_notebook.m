%MIO_NOTEBOOK  Setup for the meshio++ Octave/MATLAB example notebooks.
%   The counterpart of example/r/mio_notebook.R. Run it once, after adding
%   the binding to the path (the first cell of every notebook):
%
%     warning('off', 'Octave:shadowed-function');
%     addpath(fullfile(pwd, '..', '..', 'bindings', 'octave', 'meshioplusplus', 'inst'));
%     mio_notebook
%
%   It compiles the MEX gateway with mio_build when it is missing, and sets
%   figure defaults for inline plots. The notebook-only display helpers live
%   in their own files next to this one, one function per file so that
%   MATLAB finds them too:
%
%     render            draw a mesh with the software rasterizer, inline
%     mio_frame_image   show a frame of mio.render with image()
%     mio_bar_chart     horizontal bar chart (file sizes, cell counts, ...)
%     mio_histogram     histogram of a numeric vector
%     hex_block         an n x n x n hexahedron block fixture
%     mio_nb_temp_path  a fresh scratch file path with a given suffix
%
%   Nothing here is part of the package API.

if exist('meshioplusplus_mex', 'file') ~= 3
    % inst/*.mex is not committed: build it once against the C library
    % (MESHIOPLUSPLUS_HOME, MESHIOPLUSPLUS_CFLAGS/LIBS or pkg-config).
    mio_build();
end

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    % Without a display, Octave's only toolkit is gnuplot; its "discouraged"
    % notice would otherwise land in the first cell that draws.
    warning('off', 'Octave:gnuplot-graphics');
end
set(0, 'defaultfigurecolor', [1 1 1]);

fprintf('meshio++ %s (Octave/MATLAB binding)\n', mio.version());
