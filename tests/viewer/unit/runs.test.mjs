/**
 * The run-history helpers: filtering, sorting (unknown values last), and the
 * hyperparameter comparison that flags exactly the rows the runs disagree on.
 *
 * Run with `node --test` via tests/viewer/unit/build-and-run.mjs.
 */
import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
    COMPARED_FIELDS,
    MAX_COMPARED,
    comparisonRows,
    filterRuns,
    sortRuns,
} from '../../../src/viewer/src/dataset/runs.ts';

const run = (id, extra = {}) => ({
    job_id: id,
    status: 'finished',
    manifest: '/srv/m.json',
    tags: [],
    fields: ['q'],
    target_fields: ['T'],
    train_split: 'train',
    valid_split: 'valid',
    epochs: 100,
    batch_size: 8,
    learning_rate: 0.001,
    seed: 0,
    processor_size: 8,
    hidden_dim: 64,
    best_valid_loss: 0.001,
    final_valid_loss: 0.002,
    duration_seconds: 50,
    ...extra,
});

test('runs filter on id, status, tags and fields', () => {
    const runs = [
        run('a-1', { tags: ['sweep'] }),
        run('b-2', { status: 'failed' }),
        run('c-3', { fields: ['pressure'] }),
    ];
    assert.deepEqual(filterRuns(runs, 'b-').map((r) => r.job_id), ['b-2']);
    assert.deepEqual(filterRuns(runs, 'failed').map((r) => r.job_id), ['b-2']);
    assert.deepEqual(filterRuns(runs, 'sweep').map((r) => r.job_id), ['a-1']);
    assert.deepEqual(filterRuns(runs, 'PRESSURE').map((r) => r.job_id), ['c-3']);
    assert.equal(filterRuns(runs, '   ').length, 3);
});

test('sorting puts unknown values last in both directions', () => {
    const runs = [
        run('a', { best_valid_loss: 0.5 }),
        run('b', { best_valid_loss: null }),
        run('c', { best_valid_loss: 0.1 }),
    ];
    assert.deepEqual(
        sortRuns(runs, 'best_valid_loss', false).map((r) => r.job_id),
        ['c', 'a', 'b'],
    );
    // reversing keeps the unknown out of the "best" end
    assert.deepEqual(
        sortRuns(runs, 'best_valid_loss', true).map((r) => r.job_id),
        ['b', 'a', 'c'],
    );
    assert.deepEqual(sortRuns(runs, 'job_id', true).map((r) => r.job_id), ['c', 'b', 'a']);
    // the input is never reordered
    assert.deepEqual(runs.map((r) => r.job_id), ['a', 'b', 'c']);
});

test('a comparison flags exactly the rows that differ', () => {
    const rows = comparisonRows([
        run('a', { hidden_dim: 64, seed: 0, best_valid_loss: 0.001 }),
        run('b', { hidden_dim: 128, seed: 0, best_valid_loss: 0.002 }),
    ]);
    const by = Object.fromEntries(rows.map((r) => [r.label, r]));
    assert.equal(by.hidden.differs, true);
    assert.deepEqual(by.hidden.values, ['64', '128']);
    assert.equal(by.seed.differs, false);
    assert.equal(by['input fields'].differs, false);
    assert.deepEqual(by['input fields'].values, ['q', 'q']);
    assert.equal(by['best valid'].differs, true);
    assert.equal(rows.length, COMPARED_FIELDS.length);
    // a single run has nothing to disagree with; a missing value reads as —
    const one = comparisonRows([run('a', { duration_seconds: null })]);
    assert.equal(one.every((r) => !r.differs), true);
    assert.equal(one.find((r) => r.label === 'duration').values[0], '—');
});

test('the comparison cap is three (the validated all-pairs colour limit)', () => {
    assert.equal(MAX_COMPARED, 3);
});
