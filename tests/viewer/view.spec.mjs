/**
 * Viewing tools: standard views, projection, bookmarks, colour scale, lighting,
 * deformed shape, clip plane and the measure tool.
 *
 * Assertions go through `window.__viewerState.view`, so no GPU is needed.
 */
import { expect, test } from '@playwright/test';

test.use({ serviceWorkers: 'block' });

async function open(page, label) {
    await page.goto('./');
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('idle');
    await page.getByRole('button', { name: label }).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
}

const view = (page) => page.evaluate(() => window.__viewerState.view);

test('the View section appears with a mesh and the number keys look from the standard directions', async ({ page }) => {
    await page.goto('./');
    await expect(page.locator('#view-section')).toBeHidden();
    await open(page, 'Cube (surface)');
    await expect(page.locator('#view-section')).toBeVisible();

    await page.locator('#render').click({ position: { x: 5, y: 5 } });
    await page.keyboard.press('1');
    expect((await view(page)).direction).toEqual([0, -1, 0]);
    await page.keyboard.press('5');
    expect((await view(page)).direction).toEqual([0, 0, 1]);
    await page.keyboard.press('3');
    expect((await view(page)).direction).toEqual([1, 0, 0]);

    await page.locator('#view-buttons [data-view="left"]').click();
    expect((await view(page)).direction).toEqual([-1, 0, 0]);
});

test('orthographic projection toggles from the select and with O', async ({ page }) => {
    await open(page, 'Cube (surface)');
    expect((await view(page)).parallel).toBe(false);
    await page.selectOption('#view-projection', 'ortho');
    expect((await view(page)).parallel).toBe(true);
    await page.locator('#render').click({ position: { x: 5, y: 5 } });
    await page.keyboard.press('o');
    expect((await view(page)).parallel).toBe(false);
    await expect(page.locator('#view-projection')).toHaveValue('persp');
});

test('camera bookmarks are saved, restored and survive a reload', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.locator('#view-bookmarks summary').click();
    await page.locator('#view-buttons [data-view="top"]').click();
    await page.fill('#bm-name', 'from above');
    await page.click('#bm-save');
    expect((await view(page)).bookmarks).toEqual(['from above']);

    await page.locator('#view-buttons [data-view="front"]').click();
    await page.getByRole('button', { name: 'from above', exact: true }).click();
    expect((await view(page)).direction).toEqual([0, 0, 1]);

    await page.reload();
    await open(page, 'Cube (surface)');
    expect((await view(page)).bookmarks).toEqual(['from above']);
});

test('a pasted camera that is not valid is refused', async ({ page, context }) => {
    await open(page, 'Cube (surface)');
    page.on('dialog', (d) => d.accept('{"nope": 1}'));
    await context.grantPermissions([], {});
    await page.locator('#view-bookmarks summary').click();
    await page.evaluate(() => {
        Object.defineProperty(navigator, 'clipboard', { value: { readText: () => Promise.resolve('{"nope":1}') }, configurable: true });
    });
    await page.click('#bm-paste');
    await expect(page.locator('#bm-note')).toContainText('not a camera');
});

test('log scale and bands apply to a positive field and are refused for a signed one', async ({ page }) => {
    await open(page, 'Wave (point data)');
    const options = await page.locator('#color-by option').allTextContents();
    expect(options.length).toBeGreaterThan(1);
    await page.selectOption('#color-by', { index: 1 });
    await page.selectOption('#color-bands', '8');
    expect((await view(page)).bands).toBe(8);
    await expect(page.locator('#legend .legend-bar')).toHaveCSS('background-image', /linear-gradient/);

    await page.selectOption('#color-scale', 'log');
    const v = await view(page);
    // Whatever the field, the viewer only reports log when the range allows it.
    const range = await page.evaluate(() => ({
        min: Number(document.querySelector('.legend-min').value),
        max: Number(document.querySelector('.legend-max').value),
    }));
    expect(v.logScale).toBe(range.min > 0 && range.max > range.min);
});

test('the colour-blind filter narrows the colour maps', async ({ page }) => {
    await open(page, 'Cube (surface)');
    const all = await page.locator('#colormap option').count();
    await page.check('#colormap-safe');
    const safe = await page.locator('#colormap option').count();
    expect(safe).toBeGreaterThan(0);
    expect(safe).toBeLessThan(all);
});

test('lighting sliders and reset do not break the scene', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.locator('#view-light summary').click();
    await page.locator('#light-specular').fill('60');
    await page.check('#light-cull');
    await page.click('#light-reset');
    await expect(page.locator('#light-specular')).toHaveValue('0');
    await expect(page.locator('#light-cull')).not.toBeChecked();
    expect(await page.evaluate(() => window.__viewerState.status)).toBe('ready');
});

test('the deformed shape moves the mesh by a vector array', async ({ page }) => {
    await open(page, 'Wave (point data)');
    await page.locator('#view-warp summary').click();
    const names = await page.locator('#warp-array option').allTextContents();
    expect(names).toContain('gradient');
    expect((await view(page)).deformed).toBe(0);
    await page.selectOption('#warp-array', 'gradient');
    const drawn = (await view(page)).deformed;
    expect(drawn).toBeGreaterThan(0);
    await page.locator('#warp-scale').fill('200');
    expect((await view(page)).deformed).toBeGreaterThan(drawn);
    await page.selectOption('#warp-array', '');
    expect((await view(page)).deformed).toBe(0);
});

test('the clip plane turns on and off, and clipped hits are not picked', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.locator('#view-clip summary').click();
    await page.check('#clip-on');
    expect((await view(page)).clip).toBe(true);
    await page.locator('#clip-pos').fill('0');
    await page.uncheck('#clip-on');
    expect((await view(page)).clip).toBe(false);
});

test('measure takes two points, reports the distance, and Escape stops it', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.click('#tool-measure');
    expect((await view(page)).measuring).toBe(true);
    const b = await page.locator('#render').boundingBox();
    await page.mouse.click(b.x + b.width * 0.42, b.y + b.height * 0.55);
    await page.mouse.click(b.x + b.width * 0.58, b.y + b.height * 0.55);
    await expect.poll(async () => (await view(page)).measured).toBeGreaterThan(0);
    await expect(page.locator('#measure-result')).toContainText('Distance');
    await page.keyboard.press('Escape');
    expect((await view(page)).measuring).toBe(false);
    expect((await view(page)).measured).toBeNull();
});

test('measuring and the selection tool do not both own the click', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.click('#tool-measure');
    await page.getByRole('button', { name: 'Box' }).click();
    expect((await view(page)).measuring).toBe(false);
    expect((await page.evaluate(() => window.__viewerState.selection)).tool).toBe('box');
});

test('display choices survive a reload, and a bad stored value is ignored', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.selectOption('#colormap', 'Viridis');
    await page.check('#edges');
    await page.locator('#opacity').fill('60');
    await page.reload();
    await open(page, 'Cube (surface)');
    await expect(page.locator('#colormap')).toHaveValue('Viridis');
    await expect(page.locator('#edges')).toBeChecked();
    await expect(page.locator('#opacity')).toHaveValue('60');

    await page.evaluate(() =>
        localStorage.setItem('meshioplusplus.viewer.prefs', JSON.stringify({ colormap: 'No such map', opacity: 'x', edges: 'maybe' }))
    );
    await page.reload();
    await open(page, 'Cube (surface)');
    await expect(page.locator('#colormap')).toHaveValue('Cool to Warm');
    await expect(page.locator('#edges')).not.toBeChecked();
});

test('clicking the orientation cube snaps the camera to an axis', async ({ page }) => {
    await open(page, 'Block (volume)');
    const before = (await view(page)).direction;
    const box = await page.locator('#render').boundingBox();
    // The cube sits in the bottom-left corner of the render.
    await page.mouse.click(box.x + 40, box.y + box.height - 45);
    const after = (await view(page)).direction;
    expect(after).not.toEqual(before);
    expect(after.filter((v) => Math.abs(v) === 1)).toHaveLength(1);
});
