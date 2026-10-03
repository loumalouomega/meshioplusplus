# Benchmark trends

<script setup>
import BenchmarkTrends from './.vitepress/components/BenchmarkTrends.vue'
</script>

The weekly benchmark workflow preserves its CSVs and run metadata on the repository's `benchmark-data` branch. This page plots the published history without downloading Actions artifacts. Timings are informational, never a performance gate; hosted runners are noisy. Operation determinism checks still fail on differing output digests.

Choose a suite, an environment and a case. Cases separate mesh sizes, operations/formats, parallel backends and thread counts; environments separate CPU models, compiler versions and dependency configurations. Even the same CPU model can have different contention, so inspect repeated runs rather than attributing one spike to a commit. Missing measurements are gaps, not zeroes. Raw CSV, metadata and workflow links accompany every run.

<BenchmarkTrends />

## Storage and publication

Each successful default-branch run writes `runs/<run-id>-<attempt>/formats.{json,csv}` and `operations.{json,csv}`. Records have a versioned schema and are immutable: identical ingestion is a no-op, conflicting ingestion is an error, and a workflow rerun gets a new attempt number. The stored CSV is the original measured bytes (checksummed via `raw_csv_sha256` and validated against the JSON rows); ingestion rejects missing sidecars, checksum mismatches and row disagreements, and export verifies record locations. Operations records carry the actual native build configuration and flags per backend. Other branches retain Actions artifacts but cannot publish history. Publication is serialized and uses a job-scoped `contents: write` token; branch protection must permit that bot to update `benchmark-data`.

The existing documentation deployment exports the data into the Pages artifact and also runs after a successful trusted default-branch benchmark. There is no competing Pages deployment, external chart service or runtime GitHub API/token requirement. Before the first run, local and deployed documentation show an empty-state message. Historical Actions artifacts are not backfilled automatically because they lack the new environment metadata.

The stdlib-only `tools/bench_trends.py` supports `capture`, `ingest` and `export`. For a local preview of published data, check out `benchmark-data` into a separate directory, run `python tools/bench_trends.py export <history-directory> doc/public/benchmark-trends`, then build the docs. Never combine different environment identifiers to manufacture a longer series.
