function pipeline_run_json(json_text)
%PIPELINE_RUN_JSON  Run a settings pipeline given as JSON text.
%   mio.pipeline_run_json(JSON_TEXT) is mio.pipeline_run_file with the
%   settings document passed as a string.

    meshioplusplus_mex('mio_pipeline_run_json', mio.internal.str(json_text, 'json_text'));
end
