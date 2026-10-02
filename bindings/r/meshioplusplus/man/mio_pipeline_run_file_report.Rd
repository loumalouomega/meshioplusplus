\name{mio_pipeline_run_file_report}
\alias{mio_pipeline_run_file_report}
\alias{mio_pipeline_run_json_report}
\alias{mio_sequence_pipeline_run_file_report}
\alias{mio_sequence_pipeline_run_json_report}
\title{Pipeline structured reports as JSON}
\usage{
mio_pipeline_run_file_report(path)
mio_pipeline_run_json_report(text)
mio_sequence_pipeline_run_file_report(path)
mio_sequence_pipeline_run_json_report(text)
}
\arguments{
\item{path}{Settings file path.}
\item{text}{Settings JSON text.}
}
\value{A JSON character scalar. Requires a JSON-enabled native build.}
\description{Execute once and return JSON with steps (op plus PascalCase counters) and warnings. Non-finite counters are null. No JSON package is required; callers can parse the string with jsonlite if desired. Existing status-only APIs remain.}
