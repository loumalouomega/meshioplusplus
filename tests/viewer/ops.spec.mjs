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

/** Open Refine, fill its predicate controls, and apply. */
async function applyRefinePredicate(page, { array, compare, value }) {
    // Open the accordion by setting `open` rather than clicking the summary:
    // clicking toggles, so a helper that runs after another op has already
    // opened it would close it again.
    await page.evaluate(() => {
        document.querySelector('#op-refine-apply')?.closest('details')?.setAttribute('open', '');
    });
    await page.locator('#op-refine-array').selectOption(array);
    await page.locator('#op-refine-compare').selectOption(compare);
    await page.locator('#op-refine-value').fill(String(value));
    await page.locator('#op-refine-apply').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
}

test('refine with a predicate refines only the matching cells', async ({ page }) => {
    // The assertions are about the WIRING, not about any sample's quality
    // distribution: a threshold nothing satisfies must leave the mesh alone, and
    // one everything satisfies must give the uniform refinement. Anything in
    // between depends on the sample and would be a flaky test.
    await openSample(page, 'Cube (surface)');
    const before = (await state(page)).numCells;

    // quality first, so the predicate has a metric to threshold: composing the
    // two is what makes adaptive refinement usable from the browser.
    // `aspect_ratio`, not `scaled_jacobian`: compute_quality reports NaN where a
    // metric does not apply (and it does not for a quadrilateral in 3-D), and a
    // non-finite value deliberately never matches a predicate.
    await applyOp(page, 'Quality', 'op-quality-apply');
    expect((await state(page)).arrays).toContain('quality:aspect_ratio');

    await applyRefinePredicate(page, {
        array: 'quality:aspect_ratio',
        compare: '<',
        value: -1,
    });
    expect((await state(page)).numCells).toBe(before);

    // Undo the no-op refine, then threshold above every possible value.
    await page.locator('#ops-undo').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
    await applyRefinePredicate(page, {
        array: 'quality:aspect_ratio',
        compare: '<',
        value: 99,
    });
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

test('a derivative attaches a field without changing the geometry', async ({ page }) => {
    await openSample(page, 'Wave (point data)');
    const before = await state(page);

    await page.getByText('Derivative', { exact: true }).click();
    await page.locator('#op-grad-array').selectOption('height');
    await page.locator('#op-grad-apply').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');

    const after = await state(page);
    // Unlike every other op, this one changes no geometry -- it only adds a
    // field, which is exactly why it is worth having in the colour-by menu.
    expect(after.numCells).toBe(before.numCells);
    expect(after.numPoints).toBe(before.numPoints);
    expect(after.arrays).toContain('height:gradient (magnitude)');
    expect(after.arrays).toContain('height:gradient[2]');
    await expect(page.locator('.op-chip')).toContainText('gradient · height');

    await page.locator('#ops-undo').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
    expect((await state(page)).arrays).toEqual(before.arrays);
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

test('a section replaces the volume with its planar cross-section and undo restores it', async ({
    page,
}) => {
    await openSample(page, 'Block (volume)');
    const before = (await state(page)).numCells;

    // Section now computes the true cross-section (slice): the volume is
    // replaced by a lower-dimensional surface of section faces, not a
    // cut-away half. The count changes and stays non-empty; undo restores.
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

test('an isosurface replaces the mesh with the level set of a point field', async ({
    page,
}) => {
    await openSample(page, 'Wave (point data)');
    const before = (await state(page)).numCells;

    // The data-driven cousin of the section: the level set of `height` at 0,
    // which on this triangle sheet is a set of contour lines.
    await page.getByText('Isosurface', { exact: true }).click();
    await page.locator('#op-iso-array').selectOption('height');
    await page.locator('#op-iso-value').fill('0');
    await page.locator('#op-iso-apply').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');

    const contoured = (await state(page)).numCells;
    expect(contoured).toBeGreaterThan(0);
    expect(contoured).not.toBe(before);
    await expect(page.locator('.op-chip')).toContainText('isosurface · height = 0');

    // Undo is exact here too -- the worker replays the shortened pipeline.
    await page.locator('#ops-undo').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
    expect((await state(page)).numCells).toBe(before);
});
