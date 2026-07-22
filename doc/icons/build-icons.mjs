#!/usr/bin/env node
/**
 * Regenerates ../../src/viewer/src/ui/icons.ts from svg-ui/*.svg (each produced by
 * `make ui` from the matching tikz-ui/*.tex — see README.md).
 *
 * Pure Node: as long as svg-ui/*.svg is up to date, and it is committed for
 * exactly this reason, no TeX install is needed.
 *
 * The post-processing differs from the sibling MDPA-Preview/CAD-Preview
 * projects' because this repository has dvisvgm rather than pdftocairo — see
 * `postProcess` for the three ways their output differs.
 */
import { readFileSync, readdirSync, writeFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const SVG_DIR = path.join(HERE, 'svg-ui');
const TEX_DIR = path.join(HERE, 'tikz-ui');
const OUT_FILE = path.join(HERE, '..', '..', 'src', 'viewer', 'src', 'ui', 'icons.ts');

/**
 * Sentinel colours from preamble.tex, and the opacity each maps to.
 *
 * Drawing in the colour we actually want is impossible — the icon has to adopt
 * the surrounding text colour — and drawing in plain black would leave the
 * post-processor guessing how dvisvgm encodes it. Sentinels make the mapping
 * an exact string match.
 */
const SENTINELS = [
    { hex: '#102030', opacity: null },
    { hex: '#405060', opacity: '.6' },
    { hex: '#708090', opacity: '.3' },
];

export function postProcess(svg) {
    let out = svg
        // dvisvgm writes a single-quoted prolog and a provenance comment.
        .replace(/^<\?xml[^>]*\?>\s*/, '')
        .replace(/<!--[^]*?-->\s*/g, '')
        .trim();

    // Drop the fixed size (dvisvgm emits pt), keeping viewBox so CSS controls
    // the rendered size. Attributes may be single- or double-quoted.
    out = out.replace(/\s(?:width|height)=(['"])[^'"]*\1/g, '');

    // CRITICAL: dvisvgm wraps its output in <g id='page1'>. These SVGs are
    // inlined into one document, so eighteen icons would mean eighteen
    // elements with id="page1" -- invalid DOM, and a getElementById collision
    // in an app whose entire test hook is getElementById.
    out = out.replace(/\sid=(['"])[^'"]*\1/g, '');

    for (const { hex, opacity } of SENTINELS) {
        for (const attr of ['stroke', 'fill']) {
            const re = new RegExp(`${attr}=(['"])${hex}\\1`, 'gi');
            out = out.replace(
                re,
                opacity
                    ? `${attr}="currentColor" ${attr}-opacity="${opacity}"`
                    : `${attr}="currentColor"`
            );
        }
    }
    return out;
}

function main() {
    const svgs = readdirSync(SVG_DIR)
        .filter((f) => f.endsWith('.svg'))
        .sort();
    if (svgs.length === 0) {
        console.error('icons: svg-ui/ is empty; run `make ui` (needs pdflatex + dvisvgm)');
        process.exit(1);
    }

    // A .tex without a .svg means someone added a source and forgot `make ui`.
    const texIds = new Set(
        readdirSync(TEX_DIR)
            .filter((f) => f.endsWith('.tex'))
            .map((f) => path.basename(f, '.tex'))
    );
    const svgIds = new Set(svgs.map((f) => path.basename(f, '.svg')));
    const missing = [...texIds].filter((id) => !svgIds.has(id));
    if (missing.length) {
        console.error(
            `icons: no SVG for ${missing.join(', ')} -- run \`make ui\` (needs pdflatex).`
        );
        process.exit(1);
    }

    const entries = svgs.map((file) => {
        const id = path.basename(file, '.svg');
        const svg = postProcess(readFileSync(path.join(SVG_DIR, file), 'utf8'));
        return { id, svg };
    });

    const body = entries
        .map(({ id, svg }) => `    ${id}: ${JSON.stringify(svg)},`)
        .join('\n');

    const source = `/**
 * GENERATED FILE — do not hand-edit.
 *
 * Regenerate with:
 *   cd doc/icons && make ts
 *
 * Sources are doc/icons/tikz-ui/*.tex; see doc/icons/README.md. Each icon is inline
 * SVG using \`currentColor\`, so it takes the colour of whatever element it
 * sits in and needs no light/dark variants.
 *
 * Because \`IconId\` is a union, \`tsc --noEmit\` checks every icon reference at
 * every call site — a typo or a removed icon is a compile error.
 */

export type IconId =
${entries.map(({ id }) => `    | '${id}'`).join('\n')};

export const ICONS: Record<IconId, string> = {
${body}
};
`;

    writeFileSync(OUT_FILE, source);
    console.log(`icons: wrote ${entries.length} icons to ${path.relative(HERE, OUT_FILE)}`);
}

if (import.meta.url === `file://${process.argv[1]}`) main();
