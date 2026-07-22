#!/usr/bin/env node
/**
 * Bundle the unit tests, then run them with `node --test`.
 *
 * `node --test` cannot import TypeScript before Node 22.6 and CI is on Node
 * 20, so the `.ts` modules under test are bundled first with esbuild (a Vite
 * dependency that is already installed — this adds no package).
 *
 * This script lives under tests/viewer/, outside src/viewer/'s own
 * node_modules tree, so a bare `import 'esbuild'` cannot resolve -- Node's
 * ESM resolver only walks up from the *importing file's own* location, and
 * tests/ is a sibling of src/, not an ancestor of src/viewer/node_modules.
 * `createRequire` rooted at src/viewer/package.json reproduces the lookup
 * `npm run test:unit` (invoked with cwd=src/viewer) would do, without
 * duplicating esbuild as a tests/-local dependency.
 */
import { execFileSync } from 'node:child_process';
import { mkdtempSync, readdirSync } from 'node:fs';
import { createRequire } from 'node:module';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const viewerDir = path.resolve(here, '../../../src/viewer');
const viewerRequire = createRequire(pathToFileURL(path.join(viewerDir, 'package.json')));
const { build } = await import(pathToFileURL(viewerRequire.resolve('esbuild')).href);
const entries = readdirSync(here)
    .filter((f) => f.endsWith('.test.mjs'))
    .map((f) => path.join(here, f));

if (entries.length === 0) {
    console.error('no *.test.mjs found in tests/viewer/unit/');
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
