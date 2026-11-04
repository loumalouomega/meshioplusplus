function report = pipeline_run_file_report(path)
%PIPELINE_RUN_FILE_REPORT  Run a settings file once and return its structured report as JSON.
%   REPORT = mio.pipeline_run_file_report(PATH) executes the pipeline once and returns a char
%   row vector of JSON with the steps (op plus PascalCase counters) and the
%   warnings; non-finite counters are null. Parse it with jsondecode if
%   desired. Requires a JSON-enabled native build (mio.pipeline_has_json).
%   The status-only entry points remain.

    report = meshioplusplus_mex('mio_pipeline_run_file_report', mio.internal.str(path, 'path'));
end
