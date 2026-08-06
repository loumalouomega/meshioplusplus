"""Catalogue the preprocessed cases into a DatasetManifest, via the real CLI.

One `dataset add` per case (carrying the generator's bump parameters as entry
metadata) and one deterministic `dataset split --assign`, so the example
exercises exactly the surface it documents; the resulting
``dataset_manifest.json`` is the single source of truth every later script
(and any hand edit) works against.
"""

import argparse
import glob
import os

import numpy as np

from meshioplusplus._cli import main as cli


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases", type=str, default="cases")
    ap.add_argument("--params", type=str, default="cases_raw/params.npy")
    ap.add_argument("--manifest", type=str, default="dataset_manifest.json")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    params = np.load(args.params, allow_pickle=True)
    paths = sorted(glob.glob(os.path.join(args.cases, "case_*.vtu")))
    if len(paths) != len(params):
        raise SystemExit(
            f"{len(paths)} preprocessed cases vs {len(params)} parameter rows -- "
            "run generate_cases.py and the preprocess pipeline first"
        )
    if os.path.exists(args.manifest):
        os.remove(args.manifest)
    for path, meta in zip(paths, params):
        cli(
            [
                "dataset",
                "add",
                args.manifest,
                path,
                "--tag",
                "synthetic",
                *[f"--meta={k}={v}" for k, v in meta.items()],
            ]
        )
    cli(
        [
            "dataset",
            "split",
            args.manifest,
            "--assign",
            "train=0.8,valid=0.1,test=0.1",
            "--seed",
            str(args.seed),
        ]
    )
    cli(["dataset", "list", args.manifest])


if __name__ == "__main__":
    main()
