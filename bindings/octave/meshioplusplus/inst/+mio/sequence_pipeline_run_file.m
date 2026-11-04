function sequence_pipeline_run_file(settings_path)
%SEQUENCE_PIPELINE_RUN_FILE  Run a sequence settings document.
%   mio.sequence_pipeline_run_file(SETTINGS_PATH) runs the settings.json at
%   SETTINGS_PATH: the pipeline schema plus Mode, Input.Pattern,
%   Input.Paths, Input.Times and Input.TimeFrom. A document using none of
%   those behaves exactly as mio.pipeline_run_file.

    meshioplusplus_mex('mio_sequence_pipeline_run_file', ...
        mio.internal.str(settings_path, 'settings_path'));
end
