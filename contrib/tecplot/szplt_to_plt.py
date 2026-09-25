# MIT License -- part of meshio++ (https://github.com/loumalouomega/meshioplusplus).
"""Re-save Tecplot SZL files (``.szplt``) as binary ``.plt``, which meshio++ reads
natively, with PyTecplot (``pip install pytecplot``; it drives a licensed
Tecplot 360 in batch mode)::

    python szplt_to_plt.py run.szplt [more.szplt ...]

Each ``<name>.szplt`` becomes ``<name>.plt`` beside it. The alternative, with
no Tecplot at all, is a meshio++ build with TecIO (``doc/formats/szplt.md``).
See ``doc/routes/tecplot_szplt.md``.
"""

from __future__ import print_function

import argparse
import os


def convert(tp, path):
    out = os.path.splitext(path)[0] + ".plt"
    tp.new_layout()
    dataset = tp.data.load_tecplot_szl(path)
    tp.data.save_tecplot_plt(out, dataset=dataset)
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("files", nargs="+")
    args = parser.parse_args(argv)

    import tecplot as tp

    for path in args.files:
        print("szplt_to_plt: wrote %s" % convert(tp, path))


if __name__ == "__main__":
    main()
