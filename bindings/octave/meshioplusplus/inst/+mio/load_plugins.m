function n = load_plugins(paths)
%LOAD_PLUGINS  Load reader/writer plugins from PATHS (a path-list string); returns the number loaded.

    n = meshioplusplus_mex('mio_load_plugins', mio.internal.str(paths, 'paths'));
end
