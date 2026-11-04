function report = pipeline_run_json_report(text)
%PIPELINE_RUN_JSON_REPORT  Run settings JSON text once and return its structured report as JSON.
%   REPORT = mio.pipeline_run_json_report(TEXT) executes the pipeline once and returns a char
%   row vector of JSON with the steps (op plus PascalCase counters) and the
%   warnings; non-finite counters are null. Parse it with jsondecode if
%   desired. Requires a JSON-enabled native build (mio.pipeline_has_json).
%   The status-only entry points remain.

    report = meshioplusplus_mex('mio_pipeline_run_json_report', mio.internal.str(text, 'text'));
end
