function pipeline_run_file(settings_path)
%PIPELINE_RUN_FILE  Run a settings.json operation pipeline.
%   mio.pipeline_run_file(SETTINGS_PATH) runs a whole settings pipeline --
%   read Input.Path, apply Operations in order, write Output.Path --
%   described by a settings.json document (PascalCase ops and keys; see
%   doc/pipeline.md). Needs a library built with the JSON parser
%   (-DMESHIOPLUSPLUS_WITH_JSON=ON, the default when the
%   src/cpp/third_party/json submodule is checked out); otherwise the raised
%   error names the flag. mio.pipeline_has_json reports which build the
%   loaded library is. See also mio.pipeline_run_json,
%   mio.pipeline_run_file_report.

    meshioplusplus_mex('mio_pipeline_run_file', mio.internal.str(settings_path, 'settings_path'));
end
