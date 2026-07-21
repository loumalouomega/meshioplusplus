/**
 * The operations panel: meshio++'s own mesh operations, run in the browser.
 *
 * The two tests that matter most are `undo is exact` and `vectors survive an
 * operation`. The first proves the replay architecture — the worker keeps the
 * original bytes and re-runs a shortened pipeline, rather than trying to
 * invert anything. The second fails under any implementation that routes a
 * mesh through the flat JS representation, which is the entire reason
 * `convertSurfaceOps` exists.
 */
import { expect, test } from '@playwright/test';

async function openSample(page, label) {
    await page.goto('./');
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('idle');
    await page.getByRole('button', { name: label }).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
}

const state = (page) => page.evaluate(() => window.__viewerState);

/** Open an accordion and press its Apply. */
async function applyOp(page, summary, id) {
    await page.getByText(summary, { exact: true }).click();
    await page.locator(`#${id}`).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
}

test('quality adds per-cell metrics to the colour-by menu', async ({ page }) => {
    await openSample(page, 'Block (volume)');
    await applyOp(page, 'Quality', 'op-quality-apply');
    expect((await state(page)).arrays).toContain('quality:scaled_jacobian');
});

test('refine multiplies the cells', async ({ page }) => {
    await openSample(page, 'Cube (surface)');
    const before = (await state(page)).numCells;
    await applyOp(page, 'Refine', 'op-refine-apply');
    // A quad refines into four; the sample is 6 quads.
    expect((await state(page)).numCells).toBe(before * 4);
});

test('smooth moves nodes without changing the counts', async ({ page }) => {
    await openSample(page, 'Wave (point data)');
    const before = await state(page);
    await applyOp(page, 'Smooth', 'op-smooth-apply');
    const after = await state(page);
    expect(after.numPoints).toBe(before.numPoints);
    expect(after.numCells).toBe(before.numCells);
});

test('undo is exact', async ({ page }) => {
    // The replay architecture in one assertion: after applying an operation
    // and undoing it, the mesh must be indistinguishable from the original.
    // Nothing here is an inverse operation -- the worker re-reads the file.
    await openSample(page, 'Block (volume)');
    const before = await state(page);

    await applyOp(page, 'Refine', 'op-refine-apply');
    expect((await state(page)).numCells).not.toBe(before.numCells);

    await page.locator('#ops-undo').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');

    const after = await state(page);
    expect(after.numPoints).toBe(before.numPoints);
    expect(after.numCells).toBe(before.numCells);
    expect(after.arrays).toEqual(before.arrays);
});

test('multi-component data survives an operation', async ({ page }) => {
    // Fails under any implementation that routes the mesh through the flat JS
    // representation on the way to the operation.
    await openSample(page, 'Wave (point data)');
    expect((await state(page)).arrays).toContain('gradient[2]');

    await applyOp(page, 'Smooth', 'op-smooth-apply');
    const after = (await state(page)).arrays;
    expect(after).toContain('gradient[2]');
    expect(after).toContain('gradient (magnitude)');
});

test('operations stack and revert clears them all', async ({ page }) => {
    await openSample(page, 'Block (volume)');
    const before = await state(page);

    await applyOp(page, 'Quality', 'op-quality-apply');
    await applyOp(page, 'Refine', 'op-refine-apply');
    await expect(page.locator('.op-chip')).toHaveCount(2);

    await page.locator('#ops-revert').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
    await expect(page.locator('.op-chip')).toHaveCount(0);
    expect((await state(page)).numCells).toBe(before.numCells);
});

test('a section cuts cells away and undo restores them', async ({ page }) => {
    await openSample(page, 'Block (volume)');
    const before = (await state(page)).numCells;

    await applyOp(page, 'Section', 'op-section-apply');
    const sectioned = (await state(page)).numCells;
    expect(sectioned).toBeGreaterThan(0);
    expect(sectioned).not.toBe(before);

    await page.locator('#ops-undo').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
    expect((await state(page)).numCells).toBe(before);
});

test('partition attaches a part id to colour by', async ({ page }) => {
    await openSample(page, 'Block (volume)');
    await applyOp(page, 'Partition', 'op-partition-apply');
    expect((await state(page)).arrays).toContain('partition:part');
});
