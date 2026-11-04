function sequence_pipeline_run_json(json_text)
%SEQUENCE_PIPELINE_RUN_JSON  Run a sequence settings document given as JSON.
%   mio.sequence_pipeline_run_json(JSON_TEXT) is
%   mio.sequence_pipeline_run_file with the document as text.

    meshioplusplus_mex('mio_sequence_pipeline_run_json', mio.internal.str(json_text, 'json_text'));
end
