"""The ``pipeline`` verb: run a settings.json operation chain
(read -> operations -> write; see ``doc/pipeline.md``).

Mirrors the native CLI's verb: ``--input``/``--output`` override the two paths
in the settings file, ``--json`` prints the machine-readable report. The run
itself is :func:`meshioplusplus.run_pipeline`, the pure-Python engine -- which
is what gives this verb the per-format Python fallbacks the native CLI's
C++-only run lacks.
"""

import json


def add_args(parser):
    parser.add_argument("settings", type=str, help="the settings.json to run")
    parser.add_argument(
        "--input", type=str, default=None, help="override the settings' Input.Path"
    )
    parser.add_argument(
        "--output", type=str, default=None, help="override the settings' Output.Path"
    )
    parser.add_argument(
        "--json",
        action="store_true",
        dest="as_json",
        help="print the run report as JSON",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="print nothing on success"
    )


def _print_counter(value):
    if isinstance(value, float) and value == int(value):
        return str(int(value))
    return str(value)


def pipeline_cmd(args):
    from .. import run_pipeline

    report = run_pipeline(args.settings, input_path=args.input, output_path=args.output)
    if args.as_json:
        print(json.dumps(report, indent=2))
    elif not args.quiet:
        for i, step in enumerate(report["steps"], start=1):
            counters = ", ".join(
                f"{key}={_print_counter(value)}"
                for key, value in step.items()
                if key != "op"
            )
            print(f"step {i}: {step['op']}" + (f" ({counters})" if counters else ""))
        for warning in report["warnings"]:
            print(f"warning: {warning}")
    return 0
