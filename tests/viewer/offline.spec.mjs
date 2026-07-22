/**
 * The wheel-bundled offline page.
 *
 * `meshioplusplus.view(backend="browser")` writes a single self-contained HTML
 * file with the mesh embedded in it. This checks that such a page actually
 * renders — over `file://`, with no server and no WebAssembly — which is the
 * one property no other test covers: the Python tests verify the page's
 * *contents*, and the smoke tests exercise the *hosted* app.
 *
 * Generate the page first (see the `viewer-web` job in ci.yml):
 *
 *     python -c "...build_page(mesh)..." > /tmp/offline.html
 *     cd src/viewer && npx playwright test ../../tests/viewer/offline.spec.mjs
 *
 * Skips when no page has been generated, so a bare `npm run test:e2e` — which
 * picks up every spec in this directory — does not fail for a missing fixture.
 */
import { existsSync } from 'node:fs';

import { expect, test } from '@playwright/test';

const PAGE = process.env.VIEWER_OFFLINE_PAGE ?? '/tmp/offline.html';

test('the offline page renders its embedded mesh with no server', async ({ page }) => {
    test.skip(!existsSync(PAGE), `no generated page at ${PAGE}`);

    const errors = [];
    page.on('pageerror', (e) => errors.push(e.message));

    await page.goto(`file://${PAGE}`);
    await expect
        .poll(() => page.evaluate(() => window.__viewerState?.status), {
            message: 'the offline page never became ready',
        })
        .toBe('ready');

    const state = await page.evaluate(() => window.__viewerState);
    // The CI fixture is the 4x4x4 hexahedral block, so this pins the same
    // boundary count the hosted viewer produces -- the Python `_renderable_
    // surface` and the WASM `convertSurface` must agree.
    expect(state.numCells).toBe(96);
    expect(state.numPoints).toBe(98);
    // The sample carries a per-cell 'layer' tag; the parent-id gather must
    // bring it onto the boundary, in both the Python and the WASM path.
    expect(state.arrays).toContain('layer');
    expect(errors).toEqual([]);
});
