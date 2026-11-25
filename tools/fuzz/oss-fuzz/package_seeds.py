"""Package project-owned regression inputs, not LFS/vendor fixtures, for OSS-Fuzz."""

import argparse
import pathlib
import zipfile


def package(source, output, generated=None):
    source, output = pathlib.Path(source), pathlib.Path(output)
    for target in sorted(output.glob("meshioplusplus_fuzz_*")):
        if target.suffix or not target.is_file():
            continue
        if target.name.startswith("meshioplusplus_fuzz_read_"):
            seed_dir = target.name[len("meshioplusplus_fuzz_read_") :]
        else:
            op = target.name[len("meshioplusplus_fuzz_") :]
            seed_dir = f"op_{op}"
        files = sorted((source / seed_dir).glob("*"))
        if generated is not None:
            files += sorted((pathlib.Path(generated) / seed_dir).glob("*"))
        files = [p for p in files if p.is_file() and 0 < p.stat().st_size <= 262144]
        if not files:
            continue
        with zipfile.ZipFile(
            output / f"{target.name}_seed_corpus.zip", "w", zipfile.ZIP_DEFLATED
        ) as archive:
            for path in files:
                data = path.read_bytes()
                if data.startswith(b"version https://git-lfs.github.com/spec/v1"):
                    raise ValueError(f"LFS pointer is not a seed: {path}")
                archive.writestr(path.name, data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=pathlib.Path)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("--generated", type=pathlib.Path)
    args = parser.parse_args()
    package(args.source, args.output, args.generated)
