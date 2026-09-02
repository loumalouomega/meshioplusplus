"""Dataset manifests: a hand-editable catalogue of solution outputs.

A training run is driven off a *directory of solution outputs* -- many cases,
each possibly a time series -- not one mesh. A :class:`DatasetManifest` is the
JSON document that catalogues such a collection: each entry records the
*source plan inputs* a :class:`~meshioplusplus.TimeSeries` needs (a pattern, a
path list, explicit times, ...), a split assignment, free-form tags, an
optional group path, and open notes/metadata.

Two rules shape everything here:

- **The manifest is the single source of truth, and hand-editing it is a
  first-class workflow.** The CLI (``meshioplusplus dataset ...``), the MCP
  tools and a plain text editor all read and write the same JSON;
  serialization is stable (sorted-by-nothing: entry order is preserved, keys
  are emitted in a fixed order) so hand edits and tool edits interleave
  without clobbering each other.
- **Entries record the plan inputs, never the resolved plan.** Caching
  resolved paths/steps/times in the file would make it a hidden cache that
  goes stale the moment a case directory changes; resolution is one
  :func:`~meshioplusplus.sequence_entries` call away and stays live.

The document is a settings-family JSON (the ``doc/pipeline.md`` /
``doc/sequences.md`` conventions): PascalCase keys, lowercase enum values,
``"Version": 1``, and strict parsing -- an unknown key anywhere is an error
naming the offender, never silently ignored. See ``doc/datasets.md``.
"""

import json
import os
import random
from dataclasses import dataclass, field

from ._sequence import TimeSeries, sequence_entries

DATASET_MANIFEST_VERSION = 1

_TIME_FROM = ("auto", "file", "filename", "index")
_TOP_KEYS = ("Version", "Name", "Description", "Metadata", "Entries")
_ENTRY_KEYS = ("Id", "Source", "Split", "Tags", "Group", "Notes", "Metadata")
_SOURCE_KEYS = ("Pattern", "Path", "Paths", "Format", "Times", "TimeFrom", "Sort")


def _err(message):
    return ValueError(f"meshio++: dataset: {message}")


def _check_keys(obj, where, allowed):
    for key in obj:
        if key not in allowed:
            raise _err(f"unknown key '{key}' in {where} (known: {', '.join(allowed)})")


def _validate_source(source, where):
    """Validate a ``Source`` object; returns a normalized plain-dict copy."""
    if not isinstance(source, dict):
        raise _err(f"{where}.Source must be an object")
    _check_keys(source, f"{where}.Source", _SOURCE_KEYS)
    kinds = [k for k in ("Pattern", "Path", "Paths") if k in source]
    if len(kinds) != 1:
        raise _err(
            f"{where}.Source needs exactly one of Pattern, Path or Paths "
            f"(got {', '.join(kinds) if kinds else 'none'})"
        )
    kind = kinds[0]
    out = {}
    if kind == "Paths":
        paths = source["Paths"]
        if not isinstance(paths, (list, tuple)) or not paths:
            raise _err(f"{where}.Source.Paths must be a non-empty array of paths")
        if not all(isinstance(p, str) and p for p in paths):
            raise _err(f"{where}.Source.Paths entries must be non-empty strings")
        out["Paths"] = list(paths)
    else:
        value = source[kind]
        if not isinstance(value, str) or not value:
            raise _err(f"{where}.Source.{kind} must be a non-empty string")
        out[kind] = value
    if "Format" in source:
        if not isinstance(source["Format"], str) or not source["Format"]:
            raise _err(f"{where}.Source.Format must be a non-empty string")
        out["Format"] = source["Format"]
    if "Times" in source:
        times = source["Times"]
        if not isinstance(times, (list, tuple)) or not all(
            isinstance(t, (int, float)) and not isinstance(t, bool) for t in times
        ):
            raise _err(f"{where}.Source.Times must be an array of numbers")
        out["Times"] = [float(t) for t in times]
    if "TimeFrom" in source:
        if source["TimeFrom"] not in _TIME_FROM:
            raise _err(
                f"{where}.Source.TimeFrom must be one of "
                f"{', '.join(repr(v) for v in _TIME_FROM)}"
            )
        out["TimeFrom"] = source["TimeFrom"]
    if "Sort" in source:
        if not isinstance(source["Sort"], bool):
            raise _err(f"{where}.Source.Sort must be true or false")
        if kind != "Paths":
            # A pattern is always sorted and a single path has nothing to
            # sort; accepting the key there would imply a choice that does
            # not exist.
            raise _err(f"{where}.Source.Sort applies only with Paths")
        out["Sort"] = source["Sort"]
    return out


def _validate_metadata(value, where):
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise _err(f"{where} must be an object")
    return dict(value)


def _resolve_source_path(path, base_dir):
    if base_dir and not os.path.isabs(path):
        return os.path.join(base_dir, path)
    return path


def portable_relpath(path, base_dir):
    """``os.path.relpath(path, base_dir)``, normalized to ``/`` for storage.

    A manifest is meant to be portable, hand-editable JSON (``doc/datasets.md``),
    so a relative source written by one platform must resolve correctly when
    read back on another -- and `_glob`'s own directory/glob split
    (``_sequence.py``) takes the rightmost of either ``/`` or ``os.sep``, so a
    stored path that already uses the *native* separator on Windows works there
    too. The one thing that must never happen is writing a bare
    ``os.path.relpath()`` result: on Windows that is backslash-separated, and a
    manifest checked out and read back on Linux/macOS would then see a literal
    backslash as part of a filename, matching nothing.
    """
    return os.path.relpath(path, base_dir).replace(os.sep, "/")


@dataclass(frozen=True)
class DatasetEntry:
    """One catalogued case: a source plan plus its curation state.

    Immutable -- curation goes through :class:`DatasetManifest`, which
    replaces entries wholesale, so no half-updated entry is ever observable.
    ``base_dir`` (the directory relative source paths resolve against; set by
    :meth:`DatasetManifest.load`) is carried for resolution only and is not
    part of the document or of equality.
    """

    id: str
    source: dict
    split: str = None
    tags: tuple = ()
    group: str = None
    notes: str = None
    metadata: dict = field(default_factory=dict)
    base_dir: str = field(default=None, compare=False)

    @classmethod
    def from_dict(cls, obj, *, where="Entries[?]", base_dir=None):
        if not isinstance(obj, dict):
            raise _err(f"{where} must be an object")
        _check_keys(obj, where, _ENTRY_KEYS)
        entry_id = obj.get("Id")
        if not isinstance(entry_id, str) or not entry_id:
            raise _err(f"{where}.Id is required and must be a non-empty string")
        if "Source" not in obj:
            raise _err(f"{where}.Source is required")
        source = _validate_source(obj["Source"], where)
        split = obj.get("Split")
        if split is not None and (not isinstance(split, str) or not split):
            raise _err(f"{where}.Split must be a non-empty string")
        tags = obj.get("Tags", [])
        if not isinstance(tags, (list, tuple)) or not all(
            isinstance(t, str) and t for t in tags
        ):
            raise _err(f"{where}.Tags must be an array of non-empty strings")
        group = obj.get("Group")
        if group is not None and (not isinstance(group, str) or not group):
            raise _err(f"{where}.Group must be a non-empty string")
        notes = obj.get("Notes")
        if notes is not None and not isinstance(notes, str):
            raise _err(f"{where}.Notes must be a string")
        metadata = _validate_metadata(obj.get("Metadata"), f"{where}.Metadata")
        return cls(
            id=entry_id,
            source=source,
            split=split,
            tags=tuple(tags),
            group=group,
            notes=notes,
            metadata=metadata,
            base_dir=base_dir,
        )

    def to_dict(self):
        """The entry as its document object -- PascalCase, optional keys
        omitted when empty so a hand-edited file stays minimal."""
        out = {"Id": self.id, "Source": dict(self.source)}
        if self.split is not None:
            out["Split"] = self.split
        if self.tags:
            out["Tags"] = list(self.tags)
        if self.group is not None:
            out["Group"] = self.group
        if self.notes is not None:
            out["Notes"] = self.notes
        if self.metadata:
            out["Metadata"] = dict(self.metadata)
        return out

    def _resolved_source(self, base_dir):
        base = self.base_dir if base_dir is None else base_dir
        if "Paths" in self.source:
            return [_resolve_source_path(p, base) for p in self.source["Paths"]]
        key = "Pattern" if "Pattern" in self.source else "Path"
        return _resolve_source_path(self.source[key], base)

    def _plan_kwargs(self):
        return {
            "file_format": self.source.get("Format"),
            "times": self.source.get("Times"),
            "time_from": self.source.get("TimeFrom", "auto"),
            "sort": self.source.get("Sort", False),
        }

    def time_series(self, *, base_dir=None, **read_kwargs):
        """The entry's :class:`~meshioplusplus.TimeSeries` -- the Source
        mapped 1:1 onto the sequence machinery; plan only, no mesh read."""
        return TimeSeries(
            self._resolved_source(base_dir), **self._plan_kwargs(), **read_kwargs
        )

    def entries(self, *, base_dir=None):
        """The resolved :func:`~meshioplusplus.sequence_entries` plan
        (``{"path", "step", "time", "time_source"}`` dicts); no mesh read."""
        return sequence_entries(self._resolved_source(base_dir), **self._plan_kwargs())


class DatasetManifest:
    """An ordered, id-keyed collection of :class:`DatasetEntry`.

    See the module docstring for the two rules (single source of truth,
    plan inputs only) and ``doc/datasets.md`` for the document schema.
    """

    def __init__(
        self, entries=(), *, name=None, description=None, metadata=None, base_dir=None
    ):
        self.name = name
        self.description = description
        self.metadata = _validate_metadata(metadata, "Metadata")
        self.base_dir = base_dir
        self._entries = []
        self._by_id = {}
        for entry in entries:
            self._append(entry)

    # ------------------------------------------------------------------ #
    # Serialization                                                      #
    # ------------------------------------------------------------------ #
    @classmethod
    def load(cls, settings):
        """Load a manifest -- a dict, the JSON text itself, or a file path
        (the ``_resolve_settings`` tri-modal rule). Loading from a *file*
        records its directory as ``base_dir``, so relative source paths
        resolve against the manifest's own location and the file stays
        portable; dict/JSON-text loads resolve against the CWD."""
        base_dir = None
        if isinstance(settings, os.PathLike):
            settings = os.fspath(settings)
        if isinstance(settings, str) and not settings.lstrip().startswith("{"):
            base_dir = os.path.dirname(os.path.abspath(settings))
            with open(settings, "r", encoding="utf-8") as handle:
                doc = json.load(handle)
        elif isinstance(settings, str):
            doc = json.loads(settings)
        elif isinstance(settings, dict):
            doc = settings
        else:
            raise TypeError(
                "meshio++: dataset: the manifest must be a dict, JSON text or a path"
            )
        if not isinstance(doc, dict):
            raise _err("the manifest document must be an object")
        _check_keys(doc, "the manifest document", _TOP_KEYS)
        version = doc.get("Version", 1)
        if not isinstance(version, int) or isinstance(version, bool) or version != 1:
            raise _err(f"unsupported Version {version!r} (this build knows 1)")
        name = doc.get("Name")
        if name is not None and not isinstance(name, str):
            raise _err("Name must be a string")
        description = doc.get("Description")
        if description is not None and not isinstance(description, str):
            raise _err("Description must be a string")
        metadata = _validate_metadata(doc.get("Metadata"), "Metadata")
        raw_entries = doc.get("Entries", [])
        if not isinstance(raw_entries, list):
            raise _err("Entries must be an array")
        entries = [
            DatasetEntry.from_dict(obj, where=f"Entries[{i}]", base_dir=base_dir)
            for i, obj in enumerate(raw_entries)
        ]
        return cls(
            entries,
            name=name,
            description=description,
            metadata=metadata,
            base_dir=base_dir,
        )

    def to_dict(self):
        """The whole document, round-trip exact against :meth:`load`."""
        out = {"Version": DATASET_MANIFEST_VERSION}
        if self.name is not None:
            out["Name"] = self.name
        if self.description is not None:
            out["Description"] = self.description
        if self.metadata:
            out["Metadata"] = dict(self.metadata)
        out["Entries"] = [entry.to_dict() for entry in self._entries]
        return out

    def save(self, path):
        """Write the document (``indent=2``, trailing newline -- a stable
        serialization, so tool edits and hand edits diff cleanly). The saved
        location becomes ``base_dir``: relative sources now resolve against
        where the file actually is."""
        path = os.fspath(path)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(self.to_dict(), handle, indent=2)
            handle.write("\n")
        self.base_dir = os.path.dirname(os.path.abspath(path))
        for i, entry in enumerate(self._entries):
            self._entries[i] = self._with(entry, base_dir=self.base_dir)
            self._by_id[entry.id] = self._entries[i]

    # ------------------------------------------------------------------ #
    # Query                                                              #
    # ------------------------------------------------------------------ #
    def __len__(self):
        return len(self._entries)

    def __iter__(self):
        return iter(list(self._entries))

    def __contains__(self, entry_id):
        return entry_id in self._by_id

    def __getitem__(self, entry_id):
        try:
            return self._by_id[entry_id]
        except KeyError:
            raise KeyError(
                f"meshio++: dataset: no entry '{entry_id}' "
                f"(known: {', '.join(sorted(self._by_id)) or 'none'})"
            ) from None

    def __repr__(self):
        return f"DatasetManifest({len(self)} entr{'y' if len(self) == 1 else 'ies'})"

    def ids(self):
        """Every entry id, in entry order."""
        return [entry.id for entry in self._entries]

    def entries(self, *, split=None, tags=None, group=None):
        """The entries matching every given filter, in entry order.

        ``split`` is an exact match (an entry with no split matches only
        ``split=None``, i.e. no filter); ``tags`` means the entry must carry
        **all** of them; ``group`` matches the group path itself or any
        descendant (``"a/b"`` matches ``"a/b"`` and ``"a/b/c"``, never
        ``"a/bc"`` -- path segments, not string prefixes).
        """
        wanted_tags = set([tags] if isinstance(tags, str) else tags or ())
        selected = []
        for entry in self._entries:
            if split is not None and entry.split != split:
                continue
            if wanted_tags and not wanted_tags.issubset(entry.tags):
                continue
            if group is not None:
                have = entry.group or ""
                if have != group and not have.startswith(group + "/"):
                    continue
            selected.append(entry)
        return selected

    def splits(self):
        """``{split: entry count}``; entries with no split count under
        ``None``."""
        counts = {}
        for entry in self._entries:
            counts[entry.split] = counts.get(entry.split, 0) + 1
        return counts

    # ------------------------------------------------------------------ #
    # Curation                                                           #
    # ------------------------------------------------------------------ #
    def add(
        self,
        source,
        *,
        id=None,
        split=None,
        tags=(),
        group=None,
        notes=None,
        metadata=None,
        validate_source=True,
    ):
        """Add an entry and return it.

        ``source`` is a Source object (``{"Pattern": ...}``), a plain string
        (a pattern when it contains ``*``/``?``, else a path), or a list of
        paths. ``id=None`` derives one from the source's stem; a derived or
        given id colliding with an existing entry is an error, never a silent
        overwrite. ``validate_source=True`` expands the plan once
        (:func:`~meshioplusplus.sequence_entries`), so a glob matching
        nothing fails here, by name, not at first training access.
        """
        if isinstance(source, dict):
            source = _validate_source(source, "add()")
        elif isinstance(source, (str, os.PathLike)):
            text = os.fspath(source)
            key = "Pattern" if ("*" in text or "?" in text) else "Path"
            source = {key: text}
        else:
            source = _validate_source(
                {"Paths": [os.fspath(p) for p in source]}, "add()"
            )
        entry_id = self._derive_id(source) if id is None else id
        if not isinstance(entry_id, str) or not entry_id:
            raise _err("the entry id must be a non-empty string")
        if entry_id in self._by_id:
            raise _err(f"an entry '{entry_id}' already exists; pass a distinct id")
        entry = DatasetEntry(
            id=entry_id,
            source=source,
            split=split,
            tags=tuple([tags] if isinstance(tags, str) else tags),
            group=group,
            notes=notes,
            metadata=_validate_metadata(metadata, "metadata"),
            base_dir=self.base_dir,
        )
        if validate_source:
            entry.entries()  # raises, naming the offender, on a bad source
        self._append(entry)
        return entry

    def remove(self, entry_id):
        """Remove the entry with this id (:class:`KeyError` if absent)."""
        entry = self[entry_id]
        self._entries.remove(entry)
        del self._by_id[entry_id]

    def set_split(self, ids, split):
        """Assign ``split`` (a string, or ``None`` to clear) to the given
        ids -- one id, an iterable, or ``"all"``."""
        if split is not None and (not isinstance(split, str) or not split):
            raise _err("the split must be a non-empty string or None")
        for entry_id in self._id_list(ids):
            self._replace(entry_id, split=split)

    def assign_splits(self, fractions, *, seed=0, by_group=False):
        """Assign every entry to a split by fractions, deterministically.

        ``fractions`` is ``{"train": 0.8, "valid": 0.1, "test": 0.1}`` and
        must sum to 1. The shuffle is ``random.Random(seed)`` over the sorted
        ids, so the same manifest and seed always produce the same
        assignment. ``by_group=True`` keeps entries sharing a ``Group``
        together (whole groups are assigned, targets counted in entries) --
        the leakage guard for cases that are variations of one another.
        """
        if not isinstance(fractions, dict) or not fractions:
            raise _err("fractions must be a non-empty {split: fraction} object")
        for split_name, frac in fractions.items():
            if not isinstance(split_name, str) or not split_name:
                raise _err("every split name must be a non-empty string")
            if not isinstance(frac, (int, float)) or isinstance(frac, bool) or frac < 0:
                raise _err(f"the fraction for '{split_name}' must be a number >= 0")
        total = sum(fractions.values())
        if abs(total - 1.0) > 1e-9:
            raise _err(f"the fractions must sum to 1 (got {total!r})")

        if by_group:
            units = {}
            for entry in self._entries:
                units.setdefault(entry.group or f"\x00{entry.id}", []).append(entry.id)
            unit_list = [units[key] for key in sorted(units)]
        else:
            unit_list = [[entry_id] for entry_id in sorted(self._by_id)]
        rng = random.Random(seed)
        rng.shuffle(unit_list)

        n = len(self._entries)
        names = list(fractions)
        # Cumulative targets; the last split takes every remaining entry, so
        # rounding never drops one.
        boundaries = []
        acc = 0.0
        for split_name in names[:-1]:
            acc += fractions[split_name]
            boundaries.append(int(round(acc * n)))
        boundaries.append(n)
        index, assigned = 0, 0
        for unit in unit_list:
            while index < len(names) - 1 and assigned >= boundaries[index]:
                index += 1
            for entry_id in unit:
                self._replace(entry_id, split=names[index])
            assigned += len(unit)

    def tag(self, ids, *, add=(), remove=()):
        """Add and/or remove tags on the given ids (one id, an iterable, or
        ``"all"``). Order of a tag list is preserved; duplicates never
        accumulate."""
        add = [add] if isinstance(add, str) else list(add)
        remove = {remove} if isinstance(remove, str) else set(remove)
        for entry_id in self._id_list(ids):
            entry = self[entry_id]
            tags = [t for t in entry.tags if t not in remove]
            tags += [t for t in add if t not in tags and t not in remove]
            self._replace(entry_id, tags=tuple(tags))

    def annotate(
        self, entry_id, *, notes=None, group=None, metadata=None, drop_metadata=()
    ):
        """Update an entry's notes/group/metadata. ``None`` leaves a field
        unchanged; ``metadata`` merges key-wise; ``drop_metadata`` removes
        keys."""
        entry = self[entry_id]
        changes = {}
        if notes is not None:
            changes["notes"] = notes
        if group is not None:
            changes["group"] = group
        if metadata is not None or drop_metadata:
            merged = dict(entry.metadata)
            merged.update(_validate_metadata(metadata, "metadata"))
            for key in (
                [drop_metadata] if isinstance(drop_metadata, str) else drop_metadata
            ):
                merged.pop(key, None)
            changes["metadata"] = merged
        if changes:
            self._replace(entry_id, **changes)

    # ------------------------------------------------------------------ #
    # Internals                                                          #
    # ------------------------------------------------------------------ #
    def _append(self, entry):
        if not isinstance(entry, DatasetEntry):
            raise TypeError("meshio++: dataset: entries must be DatasetEntry objects")
        if entry.id in self._by_id:
            raise _err(f"duplicate entry id '{entry.id}'")
        self._entries.append(entry)
        self._by_id[entry.id] = entry

    @staticmethod
    def _with(entry, **changes):
        fields = {
            "id": entry.id,
            "source": entry.source,
            "split": entry.split,
            "tags": entry.tags,
            "group": entry.group,
            "notes": entry.notes,
            "metadata": entry.metadata,
            "base_dir": entry.base_dir,
        }
        fields.update(changes)
        return DatasetEntry(**fields)

    def _replace(self, entry_id, **changes):
        entry = self[entry_id]
        replacement = self._with(entry, **changes)
        self._entries[self._entries.index(entry)] = replacement
        self._by_id[entry_id] = replacement

    def _id_list(self, ids):
        if ids == "all":
            return self.ids()
        ids = [ids] if isinstance(ids, str) else list(ids)
        if not ids:
            raise _err("no entry ids given")
        for entry_id in ids:
            self[entry_id]  # KeyError, naming the known ids, if absent
        return ids

    def _derive_id(self, source):
        if "Paths" in source:
            text = source["Paths"][0]
        else:
            text = source.get("Pattern", source.get("Path"))
        stem = os.path.splitext(os.path.basename(text))[0]
        stem = stem.replace("*", "").replace("?", "").rstrip("-_.")
        if not stem:
            raise _err(f"cannot derive an id from {text!r}; pass one explicitly")
        if stem in self._by_id:
            raise _err(f"the derived id '{stem}' already exists; pass a distinct id")
        return stem
