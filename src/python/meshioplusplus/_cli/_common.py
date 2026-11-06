"""Argument blocks shared by the CLI verbs.

Most verbs read one mesh and many write one: ``INFILE [OUTFILE]`` and the
``--input-format``/``--output-format`` overrides are declared here once, in the
order every verb lists them, so a verb's ``--help`` is unchanged by using it.
"""

from .._helpers import _reader_names, _writer_names


def add_io_args(
    parser,
    out_help=None,
    in_help="mesh file to be read from",
    output_format=None,
):
    """Add ``infile``, ``outfile`` (when ``out_help`` is given) and the format flags.

    ``--output-format`` is added when there is an ``outfile``, unless
    ``output_format`` says otherwise.
    """
    if output_format is None:
        output_format = out_help is not None
    parser.add_argument("infile", type=str, help=in_help)
    if out_help is not None:
        parser.add_argument("outfile", type=str, help=out_help)
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=_reader_names(),
        help="input file format",
        default=None,
    )
    if output_format:
        parser.add_argument(
            "--output-format",
            "-o",
            type=str,
            choices=_writer_names(),
            help="output file format",
            default=None,
        )
