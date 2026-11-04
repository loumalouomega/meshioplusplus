function tf = pipeline_has_json()
%PIPELINE_HAS_JSON  True when the loaded library has the settings-pipeline JSON parser.

    tf = meshioplusplus_mex('mio_pipeline_has_json');
end
