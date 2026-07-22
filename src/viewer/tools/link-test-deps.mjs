#!/usr/bin/env node
/**
 * Symlinks tests/viewer/node_modules -> src/viewer/node_modules.
 *
 * tests/viewer/*.spec.mjs (Playwright) and tests/viewer/unit/*.test.mjs live
 * outside src/viewer/'s own npm project, so a bare `import '@playwright/test'`
 * in those files cannot resolve on its own -- Node's ESM resolver only walks
 * up from the *importing file's* location, and tests/ is a sibling of src/,
 * not an ancestor of src/viewer/node_modules. Rather than installing a second
 * copy of Playwright's devDependencies under tests/viewer/, this symlink lets
 * those files resolve into the one real install. Run automatically via
 * `postinstall` so it exists before any test script (`npm run test:unit` uses
 * a `createRequire` instead, since bundling with esbuild only needs one
 * resolved entry point rather than transparent resolution for every spec
 * file Playwright loads).
 */
import { existsSync, symlinkSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const viewerDir = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const link = path.resolve(viewerDir, '../../tests/viewer/node_modules');
const target = path.join(viewerDir, 'node_modules');

if (!existsSync(link)) {
    symlinkSync(target, link, 'dir');
}
