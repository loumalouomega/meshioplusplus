/**
 * The TS manifest model must be a faithful twin of `_dataset.py`: strict
 * parsing that names the offender, and serialization byte-parity with
 * `DatasetManifest.save` — pinned here against a fixture string produced by
 * the real Python `save()` (regenerate with the snippet in the comment
 * below if the Python serialization ever legitimately changes).
 *
 * Run with `node --test` via tests/viewer/unit/build-and-run.mjs.
 */
import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
    addEntry,
    annotateEntry,
    deriveId,
    emptyManifest,
    filterEntries,
    getEntry,
    parseManifest,
    removeEntry,
    setSplit,
    splitCounts,
    stringifyManifest,
    tagEntries,
} from '../../../src/viewer/src/dataset/manifest.ts';

/* Produced by (note the em dash — it pins the \uXXXX escaping):
 *   m = DatasetManifest(name="campaign", description="Re sweep",
 *                       metadata={"solver": "kratos"})
 *   m.add({"Pattern": "runs/c42/out_*.vtu", "Format": "vtu",
 *          "Times": [0.25, 0.5], "TimeFrom": "filename"}, id="c42",
 *         split="train", tags=["re100", "coarse"], group="cyl/laminar",
 *         notes="restarted at t=0.3 — twice",
 *         metadata={"Re": 100, "Ma": 0.3}, validate_source=False)
 *   m.add({"Paths": ["a.vtu", "b.vtu"], "Sort": True}, id="pair",
 *         validate_source=False)
 *   m.add("single.vtu", validate_source=False)
 *   m.add({"Pattern": "coarse/*.vtu"}, target={"Pattern": "fine/*.vtu"},
 *         id="sr", split="train", validate_source=False)
 *   m.save(path)
 */
const PYTHON_FIXTURE = `{
  "Version": 1,
  "Name": "campaign",
  "Description": "Re sweep",
  "Metadata": {
    "solver": "kratos"
  },
  "Entries": [
    {
      "Id": "c42",
      "Source": {
        "Pattern": "runs/c42/out_*.vtu",
        "Format": "vtu",
        "Times": [
          0.25,
          0.5
        ],
        "TimeFrom": "filename"
      },
      "Split": "train",
      "Tags": [
        "re100",
        "coarse"
      ],
      "Group": "cyl/laminar",
      "Notes": "restarted at t=0.3 \\u2014 twice",
      "Metadata": {
        "Re": 100,
        "Ma": 0.3
      }
    },
    {
      "Id": "pair",
      "Source": {
        "Paths": [
          "a.vtu",
          "b.vtu"
        ],
        "Sort": true
      }
    },
    {
      "Id": "single",
      "Source": {
        "Path": "single.vtu"
      }
    },
    {
      "Id": "sr",
      "Source": {
        "Pattern": "coarse/*.vtu"
      },
      "Target": {
        "Pattern": "fine/*.vtu"
      },
      "Split": "train"
    }
  ]
}
`;

test('round-trips a Python-saved manifest byte-for-byte', () => {
    const manifest = parseManifest(PYTHON_FIXTURE);
    assert.equal(stringifyManifest(manifest), PYTHON_FIXTURE);
});

test('the documented float gap: a Python 0.0 normalizes to 0 once', () => {
    const manifest = parseManifest(
        '{"Version": 1, "Entries": [{"Id": "x", "Source": {"Path": "p", "Times": [0.0]}}]}',
    );
    const first = stringifyManifest(manifest);
    assert.ok(first.includes('\n          0\n'));
    // ... and is stable from then on.
    assert.equal(stringifyManifest(parseManifest(first)), first);
});

test('strict parsing names the offender', () => {
    const cases = [
        [{ Version: 1, Bogus: 1, Entries: [] }, "unknown key 'Bogus'"],
        [{ Version: 2, Entries: [] }, 'unsupported Version 2'],
        [{ Entries: [{ Source: { Path: 'p' } }] }, 'Entries[0].Id'],
        [{ Entries: [{ Id: 'a', Source: { Path: 'p' }, Zz: 1 }] }, "unknown key 'Zz'"],
        [{ Entries: [{ Id: 'a', Source: {} }] }, 'exactly one of'],
        [{ Entries: [{ Id: 'a', Source: { Path: 'p', Pattern: 'q' } }] }, 'exactly one of'],
        [{ Entries: [{ Id: 'a', Source: { Path: 'p', Nope: 1 } }] }, "unknown key 'Nope'"],
        [
            { Entries: [{ Id: 'a', Source: { Path: 'p', TimeFrom: 'bogus' } }] },
            'TimeFrom',
        ],
        [
            { Entries: [{ Id: 'a', Source: { Path: 'p', Sort: true } }] },
            'Sort applies only with Paths',
        ],
        [{ Entries: [{ Id: 'a', Source: { Path: 'p' }, Tags: [1] }] }, 'Tags'],
        [
            {
                Entries: [
                    { Id: 'a', Source: { Path: 'p' } },
                    { Id: 'a', Source: { Path: 'q' } },
                ],
            },
            "duplicate entry id 'a'",
        ],
    ];
    for (const [doc, needle] of cases) {
        assert.throws(
            () => parseManifest(doc),
            (e) => e.message.startsWith('meshio++: dataset:') && e.message.includes(needle),
            `expected '${needle}' for ${JSON.stringify(doc)}`,
        );
    }
});

test('empty optional keys are omitted on save', () => {
    const manifest = emptyManifest();
    addEntry(manifest, { Path: 'a.vtu' });
    const doc = JSON.parse(stringifyManifest(manifest));
    assert.deepEqual(Object.keys(doc), ['Version', 'Entries']);
    assert.deepEqual(Object.keys(doc.Entries[0]), ['Id', 'Source']);
});

test('id derivation mirrors _derive_id', () => {
    assert.equal(deriveId({ Pattern: 'runs/out_*.vtu' }, []), 'out');
    assert.equal(deriveId({ Path: 'single.vtu' }, []), 'single');
    assert.equal(deriveId({ Paths: ['p0.vtu', 'p1.vtu'] }, []), 'p0');
    assert.throws(() => deriveId({ Pattern: '*.vtu' }, []), /cannot derive an id/);
    assert.throws(() => deriveId({ Path: 'a.vtu' }, ['a']), /already exists/);
});

test('curation: split, tags, annotate, remove', () => {
    const manifest = emptyManifest();
    addEntry(manifest, { Path: 'a.vtu' }, { id: 'a', tags: ['raw'] });
    addEntry(manifest, { Path: 'b.vtu' }, { id: 'b' });
    setSplit(manifest, ['a'], 'train');
    setSplit(manifest, 'all', null);
    setSplit(manifest, ['b'], 'test');
    assert.equal(getEntry(manifest, 'a').split, null);
    assert.equal(getEntry(manifest, 'b').split, 'test');
    assert.throws(() => setSplit(manifest, ['nope'], 'x'), /no entry 'nope'/);

    tagEntries(manifest, 'all', { add: ['v1', 'raw'] });
    assert.deepEqual(getEntry(manifest, 'a').tags, ['raw', 'v1']); // order kept, no dupes
    tagEntries(manifest, ['a'], { remove: ['raw'] });
    assert.deepEqual(getEntry(manifest, 'a').tags, ['v1']);

    annotateEntry(manifest, 'a', { notes: 'odd', metadata: { Re: 100 } });
    annotateEntry(manifest, 'a', { metadata: { Ma: 0.1 }, dropMetadata: ['Re'] });
    const a = getEntry(manifest, 'a');
    assert.equal(a.notes, 'odd');
    assert.deepEqual(a.metadata, { Ma: 0.1 });

    removeEntry(manifest, 'b');
    assert.deepEqual(
        manifest.entries.map((e) => e.id),
        ['a'],
    );
});

test('filters: exact split, all-of tags, group path segments', () => {
    const manifest = emptyManifest();
    addEntry(manifest, { Path: 'a' }, { id: 'a', split: 'train', tags: ['x', 'y'], group: 'g/h' });
    addEntry(manifest, { Path: 'b' }, { id: 'b', split: 'test', tags: ['x'], group: 'g' });
    addEntry(manifest, { Path: 'c' }, { id: 'c', group: 'gh' });
    const ids = (f) => filterEntries(manifest, f).map((e) => e.id);
    assert.deepEqual(ids({}), ['a', 'b', 'c']);
    assert.deepEqual(ids({ split: 'train' }), ['a']);
    assert.deepEqual(ids({ tags: ['x'] }), ['a', 'b']);
    assert.deepEqual(ids({ tags: ['x', 'y'] }), ['a']);
    assert.deepEqual(ids({ group: 'g' }), ['a', 'b']); // 'gh' is not a descendant
    assert.deepEqual(ids({ group: 'g/h' }), ['a']);
    assert.deepEqual(
        [...splitCounts(manifest).entries()],
        [
            ['train', 1],
            ['test', 1],
            ['', 1],
        ],
    );
});
