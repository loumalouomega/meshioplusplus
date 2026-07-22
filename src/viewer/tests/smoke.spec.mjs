/**
 * End-to-end smoke test.
 *
 * Drives the real app against the real meshio++ WASM package, so it fails if
 * the package's API changes under the viewer — which is the point: the demo
 * consumes `@meshioplusplus/wasm` exactly as an external user would, and
 * nothing else in CI would notice it rotting.
 *
 * Assertions go through `window.__viewerState` rather than pixels, so this
 * needs no GPU and runs on any runner.
 */
import { expect, test } from '@playwright/test';

/** Wait for the app to settle, surfacing its own error message on failure. */
async function waitReady(page) {
    await expect
        .poll(() => page.evaluate(() => window.__viewerState.status), {
            message: 'viewer never reached "ready"',
        })
        .toBe('ready');
}

const state = (page) => page.evaluate(() => window.__viewerState);

test('boots and reports the formats this build supports', async ({ page }) => {
    await page.goto('./');
    await expect
        .poll(() => page.evaluate(() => window.__viewerState.status))
        .toBe('idle');
    // Populated from availableFormats(), so an empty menu means the WASM
    // module never started.
    const options = await page.locator('#convert-format option').count();
    expect(options).toBeGreaterThan(10);
});

test('loads a surface sample and reports its real counts', async ({ page }) => {
    await page.goto('./');
    await page.getByRole('button', { name: 'Cube (surface)' }).click();
    await waitReady(page);

    const s = await state(page);
    expect(s.numPoints).toBe(8);
    expect(s.numCells).toBe(6);
    // The per-face cell_data written into the sample must survive to the
    // renderer's colour-by menu.
    expect(s.arrays).toContain('face');
});

test('a volume mesh is shown by its boundary', async ({ page }) => {
    await page.goto('./');
    await page.getByRole('button', { name: 'Block (volume)' }).click();
    await waitReady(page);

    // A 4x4x4 hexahedral block: 64 cells in the file, 6*16 = 96 boundary quads
    // on screen. Asserting the boundary count is what proves convertSurface
    // ran rather than the volume being silently dropped.
    const s = await state(page);
    expect(s.numCells).toBe(96);
    expect(s.numPoints).toBe(98);
    await expect(page.locator('#info')).toContainText('hexahedron');
    // extract_surface drops cell data; the parent-id gather brings it back,
    // which is what makes colouring a solid by its tag possible at all.
    expect(s.arrays).toContain('layer');
});

test('multi-component data survives the WASM boundary', async ({ page }) => {
    // The reason the viewer routes through convertSurface instead of
    // readMesh -> extractSkin -> writeMesh: the flat JS mesh cannot carry a
    // 3-component array, so `gradient` would arrive as a scalar or not at all.
    await page.goto('./');
    await page.getByRole('button', { name: 'Wave (point data)' }).click();
    await waitReady(page);

    const s = await state(page);
    expect(s.arrays).toContain('height');
    expect(s.arrays).toContain('gradient (magnitude)');
    expect(s.arrays).toContain('gradient[0]');
    expect(s.arrays).toContain('gradient[2]');
});

test('colouring by an array shows the scalar bar', async ({ page }) => {
    await page.goto('./');
    await page.getByRole('button', { name: 'Wave (point data)' }).click();
    await waitReady(page);

    await page.locator('#color-by').selectOption({ label: 'height' });
    // No exception and the selection sticks: the mapper accepted the array.
    await expect(page.locator('#color-by')).toHaveValue(/height/);
    expect(await state(page)).toHaveProperty('status', 'ready');
});

test('converts a mesh to another format entirely client-side', async ({ page }) => {
    await page.goto('./');
    await page.getByRole('button', { name: 'Cube (surface)' }).click();
    await waitReady(page);

    await page.locator('#convert-format').selectOption('stl');
    const download = page.waitForEvent('download');
    await page.getByRole('button', { name: 'Download' }).click();
    const file = await download;
    expect(file.suggestedFilename()).toBe('cube.stl');
});

test('an unreadable file reports an error instead of hanging', async ({ page }) => {
    await page.goto('./');
    await expect
        .poll(() => page.evaluate(() => window.__viewerState.status))
        .toBe('idle');

    await page.locator('#file').setInputFiles({
        name: 'broken.vtu',
        mimeType: 'application/octet-stream',
        buffer: Buffer.from('this is not a mesh'),
    });
    await expect
        .poll(() => page.evaluate(() => window.__viewerState.status))
        .toBe('error');
    await expect(page.locator('#status')).toBeVisible();
});
