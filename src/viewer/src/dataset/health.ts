/**
 * Per-entry scans and the manifest-level health summary built from them.
 *
 * This is the single home of the load-bearing rule that a `quality:*`
 * array's NaN means "metric N/A for this cell type" BY DESIGN
 * (`compute_quality`'s own convention), so quality rows are excluded from
 * every NaN/Inf-based bad-case count — here, in the summary table and in the
 * entry badges — or every entry would badge spuriously. Pure data: no DOM,
 * no worker, unit-tested in `tests/viewer/unit/health.test.mjs`.
 */

import type { ArraySummary } from '../worker/protocol';
import { splitCounts, type Manifest } from './manifest';
import type { EntryScan, ManifestHealth, SplitBalance } from './types';

/** Is a summary row a `quality:*` metric? */
export function isQualityRow(s: ArraySummary): boolean {
    return s.name.startsWith('quality:');
}

/** The quality-specific bad rules: inverted/degenerate cells present, or a
 * negative worst scaled Jacobian. */
export function qualityRowIsBad(s: ArraySummary): boolean {
    if (s.name === 'quality:inverted' || s.name === 'quality:degenerate') {
        return s.mean * s.numValues > 0;
    }
    return s.name === 'quality:scaled_jacobian' && s.min < 0;
}

/** The scan recorded for an entry that could not be staged or summarized. */
export function emptyScan(): EntryScan {
    return {
        steps: 0,
        numNan: 0,
        numInf: 0,
        numInverted: 0,
        numDegenerate: 0,
        minScaledJacobian: null,
        arrays: [],
    };
}

function countFlagged(summaries: ArraySummary[], name: string): number {
    const row = summaries.find((s) => s.name === name);
    return row ? Math.round(row.mean * row.numValues) : 0;
}

/** One entry's scan from its resolved plan length and step-0 summaries. */
export function entryScanFromSummaries(steps: number, summaries: ArraySummary[]): EntryScan {
    // The NaN/Inf lane counts DATA arrays only — a quality metric's NaN is
    // "N/A for this cell type" by design, not a bad value.
    const dataRows = summaries.filter((s) => !isQualityRow(s));
    const jacobian = summaries.find((s) => s.name === 'quality:scaled_jacobian');
    return {
        steps,
        numNan: dataRows.reduce((n, s) => n + s.numNan, 0),
        numInf: dataRows.reduce((n, s) => n + s.numInf, 0),
        numInverted: countFlagged(summaries, 'quality:inverted'),
        numDegenerate: countFlagged(summaries, 'quality:degenerate'),
        minScaledJacobian:
            jacobian && Number.isFinite(jacobian.min) ? jacobian.min : null,
        arrays: dataRows.map((s) => `${s.location}:${s.name}`).sort(),
    };
}

/** Any bad signal on one entry (an unreadable entry counts as bad). */
export function scanIsBad(scan: EntryScan): boolean {
    return (
        scan.steps === 0 ||
        scan.numNan + scan.numInf > 0 ||
        scan.numInverted > 0 ||
        scan.numDegenerate > 0 ||
        (scan.minScaledJacobian !== null && scan.minScaledJacobian < 0)
    );
}

/** `{split: count}` -> ordered balance rows, unassigned (`''`) last. */
export function splitBalanceOf(counts: Map<string, number>, total: number): SplitBalance[] {
    const rows = [...counts.entries()].map(([split, count]) => ({
        split,
        count,
        fraction: total ? count / total : 0,
    }));
    rows.sort((x, y) => {
        if (x.split === '') return 1;
        if (y.split === '') return -1;
        return x.split.localeCompare(y.split);
    });
    return rows;
}

/**
 * Aggregate the scans of one manifest. Entries without a scan contribute
 * nothing (they are simply not `scanned`); an unreadable entry
 * (`steps === 0`) is bad but takes no part in the fields-missing union, since
 * it reports no arrays at all.
 */
export function healthFromScans(
    manifest: Manifest,
    scans: Map<string, EntryScan> | Record<string, EntryScan>,
): ManifestHealth {
    const lookup = scans instanceof Map ? scans : new Map(Object.entries(scans));
    const total = manifest.entries.length;
    const health: ManifestHealth = {
        producer: 'browser',
        scanned: 0,
        total,
        numNan: 0,
        numInf: 0,
        numInverted: 0,
        numDegenerate: 0,
        minScaledJacobian: null,
        splitBalance: splitBalanceOf(splitCounts(manifest), total),
        fieldsMissing: {},
        badEntries: [],
    };
    const union = new Set<string>();
    const readable: { id: string; scan: EntryScan }[] = [];
    for (const entry of manifest.entries) {
        const scan = lookup.get(entry.id);
        if (!scan) continue;
        health.scanned += 1;
        health.numNan += scan.numNan;
        health.numInf += scan.numInf;
        health.numInverted += scan.numInverted;
        health.numDegenerate += scan.numDegenerate;
        if (scan.minScaledJacobian !== null) {
            health.minScaledJacobian =
                health.minScaledJacobian === null
                    ? scan.minScaledJacobian
                    : Math.min(health.minScaledJacobian, scan.minScaledJacobian);
        }
        if (scanIsBad(scan)) health.badEntries.push(entry.id);
        if (scan.steps > 0) {
            readable.push({ id: entry.id, scan });
            for (const a of scan.arrays) union.add(a);
        }
    }
    for (const { id, scan } of readable) {
        const have = new Set(scan.arrays);
        const missing = [...union].filter((a) => !have.has(a)).sort();
        if (missing.length) health.fieldsMissing[id] = missing;
    }
    return health;
}

export interface HealthBadge {
    label: string;
    level: 'ok' | 'warn' | 'bad';
}

/** The badges a card / the drill-down shows for a health summary. */
export function healthBadges(health: ManifestHealth): HealthBadge[] {
    const badges: HealthBadge[] = [];
    if (health.scanned < health.total) {
        badges.push({ label: `scanned ${health.scanned}/${health.total}`, level: 'warn' });
    } else if (health.total > 0) {
        badges.push({ label: `scanned ${health.total}`, level: 'ok' });
    }
    const unassigned = health.splitBalance.find((s) => s.split === '');
    if (unassigned) {
        badges.push({ label: `unassigned ${unassigned.count}`, level: 'warn' });
    }
    if (health.numNan + health.numInf > 0) {
        badges.push({ label: `NaN/Inf ${health.numNan + health.numInf}`, level: 'bad' });
    }
    if (health.numInverted > 0) {
        badges.push({ label: `inverted ${health.numInverted}`, level: 'bad' });
    }
    if (health.numDegenerate > 0) {
        badges.push({ label: `degenerate ${health.numDegenerate}`, level: 'bad' });
    }
    if (health.minScaledJacobian !== null) {
        const v = health.minScaledJacobian;
        badges.push({
            label: `min SJ ${Number(v.toPrecision(3))}`,
            level: v < 0 ? 'bad' : v < 0.2 ? 'warn' : 'ok',
        });
    }
    const missing = Object.keys(health.fieldsMissing).length;
    if (missing) {
        badges.push({ label: `fields missing in ${missing}`, level: 'warn' });
    }
    const unreadable = health.badEntries.length;
    if (unreadable && !badges.some((b) => b.level === 'bad')) {
        badges.push({ label: `unreadable ${unreadable}`, level: 'bad' });
    }
    if (!badges.some((b) => b.level !== 'ok') && health.scanned === health.total) {
        badges.push({ label: 'healthy', level: 'ok' });
    }
    return badges;
}

/** The server producer's report (`dataset_health`, snake_case) — the
 * superset shape `mcp/_health.py` emits. */
export interface ServerEntryScan {
    steps: number;
    num_nan: number;
    num_inf: number;
    num_inverted: number;
    num_degenerate: number;
    min_scaled_jacobian: number | null;
    arrays: string[];
    error?: string;
}

export interface ServerHealthReport {
    producer: 'server';
    name: string | null;
    num_entries: number;
    scanned: number;
    splits: Record<string, number>;
    split_balance: SplitBalance[];
    entries: Record<string, ServerEntryScan>;
    fields_missing: Record<string, string[]>;
    totals: {
        num_nan: number;
        num_inf: number;
        num_inverted: number;
        num_degenerate: number;
        min_scaled_jacobian: number | null;
    };
    bad_entries: string[];
    manifest_path: string;
    sha256: string;
}

/** Adapt the server's report to the browser's shapes, so the card renderer
 * and the entry badges take either producer. */
export function fromServerHealth(report: ServerHealthReport): {
    health: ManifestHealth;
    scans: Record<string, EntryScan>;
} {
    const scans: Record<string, EntryScan> = {};
    for (const [id, s] of Object.entries(report.entries)) {
        scans[id] = {
            steps: s.steps,
            numNan: s.num_nan,
            numInf: s.num_inf,
            numInverted: s.num_inverted,
            numDegenerate: s.num_degenerate,
            minScaledJacobian: s.min_scaled_jacobian,
            arrays: [...s.arrays],
        };
    }
    return {
        health: {
            producer: 'server',
            scanned: report.scanned,
            total: report.num_entries,
            numNan: report.totals.num_nan,
            numInf: report.totals.num_inf,
            numInverted: report.totals.num_inverted,
            numDegenerate: report.totals.num_degenerate,
            minScaledJacobian: report.totals.min_scaled_jacobian,
            splitBalance: report.split_balance.map((r) => ({ ...r })),
            fieldsMissing: { ...report.fields_missing },
            badEntries: [...report.bad_entries],
        },
        scans,
    };
}
