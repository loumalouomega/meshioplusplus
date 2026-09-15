/**
 * The per-entry scan and manifest-level health summary — including the
 * load-bearing rule that a `quality:*` NaN is "metric N/A", never a bad
 * value, so it stays out of every NaN/Inf count.
 *
 * Run with `node --test` via tests/viewer/unit/build-and-run.mjs.
 */
import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
    emptyScan,
    entryScanFromSummaries,
    healthBadges,
    healthFromScans,
    isQualityRow,
    qualityRowIsBad,
    scanIsBad,
    splitBalanceOf,
} from '../../../src/viewer/src/dataset/health.ts';
import { parseManifest, splitCounts } from '../../../src/viewer/src/dataset/manifest.ts';

const row = (name, extra = {}) => ({
    location: 'point_data',
    name,
    dtype: 'float64',
    numValues: 8,
    numComponents: 1,
    min: 0,
    max: 1,
    mean: 0.5,
    numNan: 0,
    numInf: 0,
    ...extra,
});

test('quality rows are recognised and judged by their own rules', () => {
    assert.equal(isQualityRow(row('quality:volume')), true);
    assert.equal(isQualityRow(row('T')), false);
    assert.equal(qualityRowIsBad(row('quality:inverted', { mean: 0.25 })), true);
    assert.equal(qualityRowIsBad(row('quality:inverted', { mean: 0 })), false);
    assert.equal(qualityRowIsBad(row('quality:scaled_jacobian', { min: -0.1 })), true);
    assert.equal(qualityRowIsBad(row('quality:scaled_jacobian', { min: 0.3 })), false);
    // a quality NaN is "N/A for this cell type", never bad
    assert.equal(qualityRowIsBad(row('quality:warpage', { numNan: 8 })), false);
});

test('a scan counts NaN/Inf over data arrays only and reads the quality lane', () => {
    const summaries = [
        row('T', { numNan: 2 }),
        row('mat', { location: 'cell_data', dtype: 'int64', numInf: 1 }),
        row('quality:warpage', { numNan: 8 }), // N/A for every cell here
        row('quality:scaled_jacobian', { min: -0.1 }),
        row('quality:inverted', { location: 'cell_data', mean: 0.25 }),
        row('quality:degenerate', { location: 'cell_data', mean: 0 }),
    ];
    const scan = entryScanFromSummaries(3, summaries);
    assert.equal(scan.steps, 3);
    assert.equal(scan.numNan, 2);
    assert.equal(scan.numInf, 1);
    assert.equal(scan.numInverted, 2);
    assert.equal(scan.numDegenerate, 0);
    assert.equal(scan.minScaledJacobian, -0.1);
    assert.deepEqual(scan.arrays, ['cell_data:mat', 'point_data:T']);
    assert.equal(scanIsBad(scan), true);
    assert.equal(scanIsBad(entryScanFromSummaries(1, [row('T')])), false);
    // an unreadable entry is bad by itself
    assert.equal(scanIsBad(emptyScan()), true);
    // a NaN-free, quality-N/A-only entry is healthy
    const clean = entryScanFromSummaries(1, [row('T'), row('quality:warpage', { numNan: 8 })]);
    assert.equal(clean.numNan, 0);
    assert.equal(clean.minScaledJacobian, null);
});

const MANIFEST = parseManifest({
    Version: 1,
    Entries: [
        { Id: 'a', Source: { Path: 'a.vtu' }, Split: 'train' },
        { Id: 'b', Source: { Path: 'b.vtu' }, Split: 'train' },
        { Id: 'c', Source: { Path: 'c.vtu' }, Split: 'valid' },
        { Id: 'd', Source: { Path: 'd.vtu' } },
    ],
});

test('split balance orders by name with unassigned last', () => {
    const rows = splitBalanceOf(splitCounts(MANIFEST), 4);
    assert.deepEqual(rows, [
        { split: 'train', count: 2, fraction: 0.5 },
        { split: 'valid', count: 1, fraction: 0.25 },
        { split: '', count: 1, fraction: 0.25 },
    ]);
});

test('the health summary aggregates scans, finds missing fields and bad entries', () => {
    const scans = new Map([
        ['a', entryScanFromSummaries(1, [row('T'), row('q'), row('quality:scaled_jacobian', { min: 0.4 })])],
        ['b', entryScanFromSummaries(2, [row('T', { numNan: 3 }), row('quality:scaled_jacobian', { min: 0.1 })])],
        ['c', emptyScan()], // unreadable
        // 'd' is not scanned
    ]);
    const health = healthFromScans(MANIFEST, scans);
    assert.equal(health.producer, 'browser');
    assert.equal(health.scanned, 3);
    assert.equal(health.total, 4);
    assert.equal(health.numNan, 3);
    assert.equal(health.minScaledJacobian, 0.1);
    // b lacks q, which a carries; the unreadable c takes no part in the union
    assert.deepEqual(health.fieldsMissing, { b: ['point_data:q'] });
    assert.deepEqual(health.badEntries, ['b', 'c']);
    assert.equal(health.splitBalance.at(-1).split, '');

    const labels = healthBadges(health).map((b) => `${b.level}:${b.label}`);
    assert.ok(labels.includes('warn:scanned 3/4'));
    assert.ok(labels.includes('warn:unassigned 1'));
    assert.ok(labels.includes('bad:NaN/Inf 3'));
    assert.ok(labels.includes('warn:min SJ 0.1'));
    assert.ok(labels.includes('warn:fields missing in 1'));
    assert.ok(!labels.some((l) => l === 'ok:healthy'));
});

test('a fully scanned clean manifest is simply healthy', () => {
    const manifest = parseManifest({
        Version: 1,
        Entries: [{ Id: 'a', Source: { Path: 'a.vtu' }, Split: 'train' }],
    });
    const scans = { a: entryScanFromSummaries(1, [row('T'), row('quality:scaled_jacobian', { min: 0.9 })]) };
    const badges = healthBadges(healthFromScans(manifest, scans));
    assert.deepEqual(
        badges.map((b) => b.level),
        ['ok', 'ok', 'ok'],
    );
    assert.equal(badges.at(-1).label, 'healthy');
});

test('the server report adapts to the browser shapes', async () => {
    const { fromServerHealth } = await import('../../../src/viewer/src/dataset/health.ts');
    const report = {
        producer: 'server',
        name: 'h',
        num_entries: 2,
        scanned: 2,
        splits: { train: 1, '': 1 },
        split_balance: [
            { split: 'train', count: 1, fraction: 0.5 },
            { split: '', count: 1, fraction: 0.5 },
        ],
        entries: {
            a: { steps: 3, num_nan: 0, num_inf: 0, num_inverted: 0, num_degenerate: 0, min_scaled_jacobian: 0.7, arrays: ['point_data:T'], target_steps: 3 },
            b: { steps: 0, num_nan: 0, num_inf: 0, num_inverted: 0, num_degenerate: 0, min_scaled_jacobian: null, arrays: [], error: 'gone' },
        },
        fields_missing: {},
        totals: { num_nan: 0, num_inf: 0, num_inverted: 0, num_degenerate: 0, min_scaled_jacobian: 0.7 },
        bad_entries: ['b'],
        manifest_path: '/w/m.json',
        sha256: 'abc',
    };
    const { health, scans } = fromServerHealth(report);
    assert.equal(health.producer, 'server');
    assert.equal(health.total, 2);
    assert.equal(health.minScaledJacobian, 0.7);
    assert.deepEqual(health.badEntries, ['b']);
    assert.deepEqual(scans.a, { steps: 3, numNan: 0, numInf: 0, numInverted: 0, numDegenerate: 0, minScaledJacobian: 0.7, arrays: ['point_data:T'], targetSteps: 3, pairingError: null });
    assert.equal(scans.b.steps, 0);
    // an entry the server said nothing about pairing for is self-supervised,
    // not broken -- both fields fall back to null rather than undefined
    assert.equal(scans.b.targetSteps, null);
    assert.equal(scans.b.pairingError, null);
});
