#!/usr/bin/env node
/**
 * CI gate for the icon set. Needs no TeX install.
 *
 * Two failure modes, both silent otherwise:
 *   1. a `.tex` was added or edited but `make ui` was never run, so the
 *      committed SVG (and therefore the generated TS) is stale;
 *   2. `build-icons.mjs` changed but `make ts` was never run, so the committed
 *      icons.ts no longer matches what the generator would produce.
 *
 * Also asserts the invariants the app depends on: `currentColor` so icons
 * theme themselves, and no `id` attribute — dvisvgm emits `<g id='page1'>`,
 * and eighteen of those inlined into one document would collide.
 */
import { readFileSync, readdirSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import { postProcess } from './build-icons.mjs';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const problems = [];

const texIds = readdirSync(path.join(HERE, 'tikz-ui'))
    .filter((f) => f.endsWith('.tex'))
    .map((f) => path.basename(f, '.tex'))
    .sort();
const svgIds = readdirSync(path.join(HERE, 'svg-ui'))
    .filter((f) => f.endsWith('.svg'))
    .map((f) => path.basename(f, '.svg'))
    .sort();

for (const id of texIds) {
    if (!svgIds.includes(id)) {
        problems.push(`tikz-ui/${id}.tex has no svg-ui/${id}.svg — run \`make ui\``);
    }
}
for (const id of svgIds) {
    if (!texIds.includes(id)) {
        problems.push(`svg-ui/${id}.svg has no tikz-ui/${id}.tex — stale output?`);
    }
}

const generated = readFileSync(
    path.join(HERE, '..', '..', 'src', 'viewer', 'src', 'ui', 'icons.ts'),
    'utf8'
);

for (const id of svgIds) {
    const svg = postProcess(readFileSync(path.join(HERE, 'svg-ui', `${id}.svg`), 'utf8'));
    if (!generated.includes(JSON.stringify(svg))) {
        problems.push(`icons.ts is stale for '${id}' — run \`cd doc/icons && make ts\``);
    }
    if (!svg.includes('currentColor')) {
        problems.push(`${id}: no currentColor — did it draw outside the sentinel colours?`);
    }
    if (/\sid=/.test(svg)) {
        problems.push(`${id}: an id attribute survived post-processing`);
    }
}

if (problems.length) {
    console.error('icon check failed:');
    for (const p of problems) console.error(`  - ${p}`);
    process.exit(1);
}
console.log(`icons: ${svgIds.length} icons are current and well-formed`);
