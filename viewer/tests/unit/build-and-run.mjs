#!/usr/bin/env node
/**
 * Bundle the unit tests, then run them with `node --test`.
 *
 * `node --test` cannot import TypeScript before Node 22.6 and CI is on Node
 * 20, so the `.ts` modules under test are bundled first with esbuild (a Vite
 * dependency that is already installed — this adds no package).
 */
import { execFileSync } from 'node:child_process';
import { mkdtempSync, readdirSync } from 'node:fs';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import { build } from 'esbuild';

const here = path.dirname(fileURLToPath(import.meta.url));
const entries = readdirSync(here)
    .filter((f) => f.endsWith('.test.mjs'))
    .map((f) => path.join(here, f));

if (entries.length === 0) {
    console.error('no *.test.mjs found in tests/unit/');
    process.exit(1);
}

const outdir = mkdtempSync(path.join(tmpdir(), 'meshio-viewer-unit-'));
await build({
    entryPoints: entries,
    outdir,
    bundle: true,
    format: 'esm',
    platform: 'node',
    target: 'node18',
    // The temp dir has no package.json, so Node would read a bare `.js` as
    // CommonJS and reject the ESM output.
    outExtension: { '.js': '.mjs' },
    // node:test and friends must stay external or esbuild inlines nothing useful.
    external: ['node:*'],
    logLevel: 'warning',
});

execFileSync(process.execPath, ['--test', outdir], { stdio: 'inherit' });
