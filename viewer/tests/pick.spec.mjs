/**
 * Click-to-inspect.
 *
 * Asserts on `window.__viewerState.picked` rather than on pixels, like every
 * other spec here, so it needs no GPU. The picking itself is a CPU ray/cell
 * intersection over the polydata, so it does not depend on the renderer at all
 * — only on the scene's camera, which is why the click has to land on the mesh.
 */
import { expect, test } from '@playwright/test';

async function openSample(page, label) {
    await page.goto('./');
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('idle');
    await page.getByRole('button', { name: label }).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
}

/** Click the middle of the render area, where the fitted mesh is. */
async function clickMesh(page) {
    const box = await page.locator('#render').boundingBox();
    await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
}

const picked = (page) => page.evaluate(() => window.__viewerState.picked);

test('inspecting is off until asked for', async ({ page }) => {
    await openSample(page, 'Cube (surface)');
    await expect(page.locator('#tool-inspect')).toHaveAttribute('aria-pressed', 'false');
    await clickMesh(page);
    expect(await picked(page)).toBeNull();
    await expect(page.locator('#pick-popover')).toBeHidden();
});

test('clicking a cell reports it and its data', async ({ page }) => {
    await openSample(page, 'Cube (surface)');
    await page.locator('#tool-inspect').click();
    await expect(page.locator('#tool-inspect')).toHaveAttribute('aria-pressed', 'true');

    await clickMesh(page);
    const info = await picked(page);
    expect(info).not.toBeNull();
    expect(info.cellId).toBeGreaterThanOrEqual(0);
    expect(info.cellId).toBeLessThan(6);
    expect(info.cellType).toBe('quad');
    // The sample carries a per-face `face` tag.
    expect(info.cellValues.map((v) => v.name)).toContain('face');
    await expect(page.locator('#pick-popover')).toBeVisible();
    await expect(page.locator('#pick-body')).toContainText('Surface cell');
});

test('the nearest point is one of the picked cell\'s own vertices', async ({ page }) => {
    await openSample(page, 'Cube (surface)');
    await page.locator('#tool-inspect').click();
    await clickMesh(page);

    const info = await picked(page);
    expect(info.pointId).toBeGreaterThanOrEqual(0);
    expect(info.pointId).toBeLessThan(8);
    expect(info.position).toHaveLength(3);
    // The cube spans the unit box, so every coordinate must be 0 or 1.
    for (const c of info.position) expect([0, 1]).toContain(c);
});

test('multi-component point data reports every component', async ({ page }) => {
    await openSample(page, 'Wave (point data)');
    await page.locator('#tool-inspect').click();
    await clickMesh(page);

    const info = await picked(page);
    expect(info).not.toBeNull();
    const gradient = info.pointValues.find((v) => v.name === 'gradient');
    expect(gradient).toBeDefined();
    // Not the magnitude -- the actual three components.
    expect(gradient.components).toHaveLength(3);
});

test('provenance names the original cell of a solid', async ({ page }) => {
    // block.vtu is 64 hexahedra shown as 96 boundary quads, so the surface cell
    // id and the originating volume cell id are genuinely different numbers.
    await openSample(page, 'Block (volume)');
    await page.locator('#tool-inspect').click();
    await clickMesh(page);

    const info = await picked(page);
    expect(info).not.toBeNull();
    expect(info.parentCell).not.toBeNull();
    expect(info.parentCell).toBeGreaterThanOrEqual(0);
    expect(info.parentCell).toBeLessThan(64);
    await expect(page.locator('#pick-body')).toContainText('Original cell');
});

test('after an operation the id is labelled as post-operation', async ({ page }) => {
    // The provenance array names a cell of the *operated* mesh, not the file,
    // so the label has to say which rather than let it read as the original.
    await openSample(page, 'Block (volume)');
    await page.getByText('Quality', { exact: true }).click();
    await page.locator('#op-quality-apply').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');

    await page.locator('#tool-inspect').click();
    await clickMesh(page);
    expect(await picked(page)).not.toBeNull();
    await expect(page.locator('#pick-body')).toContainText('after ops');
});

test('loading another mesh clears the previous pick', async ({ page }) => {
    await openSample(page, 'Cube (surface)');
    await page.locator('#tool-inspect').click();
    await clickMesh(page);
    expect(await picked(page)).not.toBeNull();

    await page.getByRole('button', { name: 'Wave (point data)' }).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
    expect(await picked(page)).toBeNull();
    await expect(page.locator('#pick-popover')).toBeHidden();
});

test('turning inspect off dismisses the popover', async ({ page }) => {
    await openSample(page, 'Cube (surface)');
    await page.locator('#tool-inspect').click();
    await clickMesh(page);
    await expect(page.locator('#pick-popover')).toBeVisible();

    await page.locator('#pick-close').click();
    await expect(page.locator('#pick-popover')).toBeHidden();
    await expect(page.locator('#tool-inspect')).toHaveAttribute('aria-pressed', 'false');
    expect(await picked(page)).toBeNull();
});
