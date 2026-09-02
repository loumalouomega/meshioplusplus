"""The nested `dataset` subcommand group: curate a dataset manifest.

`meshioplusplus dataset add / list / split / tag / annotate` over the
hand-editable `DatasetManifest` JSON (`doc/datasets.md`). The `data` group's
wiring pattern verbatim: each inner parser's `set_defaults(func=...)`
overrides the outer one, so `_main.py`'s single `args.func(args)` dispatch
needs no change. Python CLI only, like `data export` — the manifest never
reaches the C++ core and has no native-CLI counterpart.

Every mutating verb is load -> mutate -> save: the manifest file is the
single source of truth, so hand edits made between two CLI calls survive.
Sources given on the command line are stored relative to the manifest's own
directory (absolute paths stay absolute), which is what keeps the manifest
portable — see `doc/datasets.md`.
"""

import json
import os

from .._dataset import DatasetManifest, portable_relpath


def _load_or_new(path):
    if os.path.exists(path):
        return DatasetManifest.load(path)
    return DatasetManifest(base_dir=os.path.dirname(os.path.abspath(path)))


def _meta_pairs(pairs):
    """`K=V` options -> a metadata dict; V parses as JSON when it can, so
    `--meta Re=100` is the number 100 and `--meta note='"x"'` a string."""
    meta = {}
    for pair in pairs or ():
        key, sep, value = pair.partition("=")
        if not sep or not key:
            raise SystemExit(f"error: --meta expects KEY=VALUE, got '{pair}'")
        try:
            meta[key] = json.loads(value)
        except json.JSONDecodeError:
            meta[key] = value
    return meta


def _manifest_relative(source, manifest_path):
    """Store a CWD-relative command-line source relative to the manifest's
    directory (the resolution anchor); absolute paths stay absolute."""
    if os.path.isabs(source):
        return source
    base = os.path.dirname(os.path.abspath(manifest_path))
    return portable_relpath(os.path.abspath(source), base)


def _split_csv(text):
    return [item for item in (text or "").split(",") if item]


# --- add -------------------------------------------------------------------
def add_add_args(parser):
    parser.add_argument("manifest", type=str, help="manifest JSON (created if absent)")
    parser.add_argument(
        "source",
        type=str,
        nargs="+",
        help="one glob/path (quote the glob), or several paths forming one case",
    )
    parser.add_argument("--id", type=str, default=None, help="entry id (default: stem)")
    parser.add_argument("--format", type=str, default=None, help="forced input format")
    parser.add_argument(
        "--time-from",
        type=str,
        choices=["auto", "file", "filename", "index"],
        default=None,
        help="time-value source (doc/sequences.md precedence)",
    )
    parser.add_argument(
        "--times", type=str, default=None, help="explicit times, comma-separated"
    )
    parser.add_argument(
        "--sort",
        action="store_true",
        help="natural-numeric sort an explicit path list",
    )
    parser.add_argument("--split", type=str, default=None, help="split assignment")
    parser.add_argument("--tag", action="append", default=[], help="tag (repeatable)")
    parser.add_argument("--group", type=str, default=None, help="group path (a/b/c)")
    parser.add_argument("--notes", type=str, default=None, help="free-text notes")
    parser.add_argument(
        "--meta",
        action="append",
        default=[],
        metavar="K=V",
        help="metadata entry (repeatable; V parses as JSON when it can)",
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="skip expanding the source now (accept a not-yet-existing case)",
    )


def add_cmd(args):
    manifest = _load_or_new(args.manifest)
    sources = [_manifest_relative(s, args.manifest) for s in args.source]
    if len(sources) == 1:
        source = sources[0]
        key = "Pattern" if ("*" in source or "?" in source) else "Path"
        source = {key: source}
    else:
        source = {"Paths": sources}
    if args.format:
        source["Format"] = args.format
    if args.times:
        source["Times"] = [float(t) for t in _split_csv(args.times)]
    if args.time_from:
        source["TimeFrom"] = args.time_from
    if args.sort:
        source["Sort"] = True
    entry = manifest.add(
        source,
        id=args.id,
        split=args.split,
        tags=args.tag,
        group=args.group,
        notes=args.notes,
        metadata=_meta_pairs(args.meta) or None,
        validate_source=not args.no_validate,
    )
    manifest.save(args.manifest)
    steps = "not validated" if args.no_validate else f"{len(entry.entries())} step(s)"
    print(f"added '{entry.id}' ({steps}); {len(manifest)} entr(ies) total")
    return 0


# --- list ------------------------------------------------------------------
def add_list_args(parser):
    parser.add_argument("manifest", type=str, help="manifest JSON to read")
    parser.add_argument("--split", type=str, default=None, help="filter by split")
    parser.add_argument(
        "--tag", action="append", default=[], help="filter: must carry ALL (repeatable)"
    )
    parser.add_argument(
        "--group", type=str, default=None, help="filter by group path (or descendant)"
    )
    parser.add_argument(
        "--resolve",
        action="store_true",
        help="expand each entry's plan (checks files exist; reads no mesh)",
    )
    parser.add_argument("--json", action="store_true", help="emit as JSON")


def list_cmd(args):
    manifest = DatasetManifest.load(args.manifest)
    entries = manifest.entries(
        split=args.split, tags=args.tag or None, group=args.group
    )

    if args.json:
        payload = []
        for entry in entries:
            item = entry.to_dict()
            if args.resolve:
                item["Resolved"] = entry.entries()
            payload.append(item)
        print(json.dumps(payload, indent=2))
        return 0

    header = f"<meshio++ dataset> ({len(entries)}"
    if len(entries) != len(manifest):
        header += f" of {len(manifest)}"
    print(header + ")")
    if manifest.name:
        print(f"  name: {manifest.name}")
    if not entries:
        print("  No matching entries.")
        return 0
    for entry in entries:
        bits = []
        if entry.split:
            bits.append(entry.split)
        if entry.tags:
            bits.append("tags: " + ",".join(entry.tags))
        if entry.group:
            bits.append(f"group: {entry.group}")
        if args.resolve:
            bits.append(f"{len(entry.entries())} step(s)")
        source = entry.source.get(
            "Pattern", entry.source.get("Path", entry.source.get("Paths"))
        )
        suffix = f" [{'; '.join(bits)}]" if bits else ""
        print(f"  {entry.id}: {source}{suffix}")
    counts = manifest.splits()
    if set(counts) != {None}:
        summary = ", ".join(
            f"{name or '(unassigned)'}={count}" for name, count in counts.items()
        )
        print(f"  splits: {summary}")
    return 0


# --- split -----------------------------------------------------------------
def add_split_args(parser):
    parser.add_argument("manifest", type=str, help="manifest JSON to update")
    parser.add_argument(
        "--id", action="append", default=[], help="entry id (repeatable)"
    )
    parser.add_argument("--all", action="store_true", help="every entry")
    parser.add_argument(
        "--set", type=str, default=None, help="assign this split to the given ids"
    )
    parser.add_argument(
        "--assign",
        type=str,
        default=None,
        metavar="S=F,S=F",
        help="assign every entry by fractions, e.g. train=0.8,valid=0.1,test=0.1",
    )
    parser.add_argument("--seed", type=int, default=0, help="shuffle seed (--assign)")
    parser.add_argument(
        "--by-group",
        action="store_true",
        help="keep entries sharing a Group together (--assign)",
    )


def split_cmd(args):
    if (args.set is None) == (args.assign is None):
        print("error: pass exactly one of --set or --assign")
        return 2
    manifest = DatasetManifest.load(args.manifest)
    if args.set is not None:
        ids = "all" if args.all else args.id
        if not ids:
            print("error: --set needs --id (repeatable) or --all")
            return 2
        manifest.set_split(ids, args.set)
    else:
        fractions = {}
        for pair in _split_csv(args.assign):
            name, sep, frac = pair.partition("=")
            if not sep:
                print(f"error: --assign expects SPLIT=FRACTION, got '{pair}'")
                return 2
            fractions[name] = float(frac)
        manifest.assign_splits(fractions, seed=args.seed, by_group=args.by_group)
    manifest.save(args.manifest)
    counts = manifest.splits()
    print(
        ", ".join(f"{name or '(unassigned)'}={count}" for name, count in counts.items())
    )
    return 0


# --- tag -------------------------------------------------------------------
def add_tag_args(parser):
    parser.add_argument("manifest", type=str, help="manifest JSON to update")
    parser.add_argument(
        "--id", action="append", default=[], help="entry id (repeatable)"
    )
    parser.add_argument("--all", action="store_true", help="every entry")
    parser.add_argument("--add", type=str, default=None, help="tags to add (comma-sep)")
    parser.add_argument(
        "--remove", type=str, default=None, help="tags to remove (comma-sep)"
    )


def tag_cmd(args):
    ids = "all" if args.all else args.id
    if not ids:
        print("error: pass --id (repeatable) or --all")
        return 2
    if args.add is None and args.remove is None:
        print("error: pass --add and/or --remove")
        return 2
    manifest = DatasetManifest.load(args.manifest)
    manifest.tag(ids, add=_split_csv(args.add), remove=_split_csv(args.remove))
    manifest.save(args.manifest)
    shown = manifest.ids() if ids == "all" else ids
    for entry_id in shown:
        print(f"  {entry_id}: {','.join(manifest[entry_id].tags) or '(none)'}")
    return 0


# --- annotate --------------------------------------------------------------
def add_annotate_args(parser):
    parser.add_argument("manifest", type=str, help="manifest JSON to update")
    parser.add_argument("--id", type=str, required=True, help="entry id")
    parser.add_argument("--notes", type=str, default=None, help="replace the notes")
    parser.add_argument("--group", type=str, default=None, help="set the group path")
    parser.add_argument(
        "--meta",
        action="append",
        default=[],
        metavar="K=V",
        help="merge a metadata entry (repeatable; V parses as JSON when it can)",
    )
    parser.add_argument(
        "--del-meta", action="append", default=[], help="drop a metadata key"
    )


def annotate_cmd(args):
    manifest = DatasetManifest.load(args.manifest)
    manifest.annotate(
        args.id,
        notes=args.notes,
        group=args.group,
        metadata=_meta_pairs(args.meta) or None,
        drop_metadata=args.del_meta,
    )
    manifest.save(args.manifest)
    entry = manifest[args.id]
    print(
        f"  {entry.id}: notes={entry.notes!r} group={entry.group!r} "
        f"metadata={entry.metadata!r}"
    )
    return 0


# --- group wiring ----------------------------------------------------------
_VERBS = (
    (
        "add",
        "Add a case (a glob, one file, or several paths) to a dataset manifest "
        "(created if absent)",
        add_add_args,
        add_cmd,
    ),
    (
        "list",
        "List a manifest's entries (filter by split/tag/group)",
        add_list_args,
        list_cmd,
    ),
    (
        "split",
        "Assign splits: --set on chosen ids, or --assign by fractions "
        "(deterministic, seeded)",
        add_split_args,
        split_cmd,
    ),
    ("tag", "Add/remove tags on chosen ids (or --all)", add_tag_args, tag_cmd),
    (
        "annotate",
        "Set an entry's notes / group / metadata",
        add_annotate_args,
        annotate_cmd,
    ),
)


def add_args(parser):
    sub = parser.add_subparsers(title="dataset subcommands", dest="dataset_command")
    sub.required = True
    for name, help_text, add, func in _VERBS:
        p = sub.add_parser(name, help=help_text)
        add(p)
        p.set_defaults(func=func)


def dataset_cmd(args):
    """Reached only if argparse somehow admits a bare ``dataset`` invocation."""
    print("error: 'dataset' requires a subcommand: " + ", ".join(v[0] for v in _VERBS))
    return 2
