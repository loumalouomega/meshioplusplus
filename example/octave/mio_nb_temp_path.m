function path = mio_nb_temp_path(suffix)
%MIO_NB_TEMP_PATH  A fresh scratch file path with the given suffix.
%   PATH = mio_nb_temp_path('.vtu') names a new file in a scratch folder
%   made once per session, for one-shot round trips. Notebook helper, not
%   part of the package API.

    persistent folder counter
    if isempty(folder)
        folder = tempname();
        mkdir(folder);
        counter = 0;
    end
    counter = counter + 1;
    path = fullfile(folder, sprintf('meshioplusplus_octave_nb_%d%s', counter, suffix));
end
