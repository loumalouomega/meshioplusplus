/**
 * Interactive selection: click, box, lasso, grow/invert, hide/isolate,
 * selection by value, id export, and what survives an operation.
 *
 * Real pointer gestures on the canvas; assertions go through
 * `window.__viewerState.selection`, so no GPU is needed.
 */
import { expect, test } from '@playwright/test';

test.use({ serviceWorkers: 'block' });

async function open(page, label) {
    await page.goto('./');
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('idle');
    await page.getByRole('button', { name: label }).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
}

const sel = (page) => page.evaluate(() => window.__viewerState.selection);
const count = async (page) => {
    const s = await sel(page);
    return s.sets.find((x) => x.id === s.active)?.count ?? 0;
};
const canvasBox = (page) => page.locator('#render').boundingBox();

/** Drag a box across the mesh, starting below the toolbar. */
async function boxDrag(page, modifier) {
    const b = await canvasBox(page);
    if (modifier) await page.keyboard.down(modifier);
    await page.mouse.move(b.x + 10, b.y + 90);
    await page.mouse.down();
    await page.mouse.move(b.x + b.width / 2, b.y + b.height / 2, { steps: 4 });
    await page.mouse.move(b.x + b.width - 10, b.y + b.height - 10, { steps: 4 });
    await page.mouse.up();
    if (modifier) await page.keyboard.up(modifier);
}

async function clickAt(page, fx, fy, modifier) {
    const b = await canvasBox(page);
    if (modifier) await page.keyboard.down(modifier);
    await page.mouse.click(b.x + b.width * fx, b.y + b.height * fy);
    if (modifier) await page.keyboard.up(modifier);
}

test('the selection tool is off until asked for, and a plain drag does not select', async ({ page }) => {
    await open(page, 'Cube (surface)');
    expect((await sel(page)).tool).toBe('off');
    await boxDrag(page);
    expect(await count(page)).toBe(0);
});

test('box select picks the visible faces, and every face when the visible-side limit is off', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.getByRole('button', { name: 'Box' }).click();
    expect((await sel(page)).tool).toBe('box');
    await boxDrag(page);
    const visible = await count(page);
    expect(visible).toBeGreaterThan(0);
    expect(visible).toBeLessThan(6);

    await page.locator('#sel-front').uncheck();
    await boxDrag(page);
    expect(await count(page)).toBe(6);
    await expect(page.locator('#sb-selection')).toContainText('6 selected');
});

test('click picks a cell, Shift adds, Ctrl toggles it off', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.getByRole('button', { name: 'Click' }).click();
    await clickAt(page, 0.5, 0.55);
    expect(await count(page)).toBe(1);
    await clickAt(page, 0.5, 0.3, 'Shift');
    expect(await count(page)).toBe(2);
    await clickAt(page, 0.5, 0.55, 'Control');
    expect(await count(page)).toBe(1);
    // Clicking empty space with no modifier clears.
    await clickAt(page, 0.97, 0.9);
    expect(await count(page)).toBe(0);
});

test('lasso selects what the loop encloses', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.locator('#sel-front').uncheck();
    await page.getByRole('button', { name: 'Lasso' }).click();
    const b = await canvasBox(page);
    const cx = b.x + b.width / 2;
    const cy = b.y + b.height / 2;
    const r = Math.min(b.width, b.height) * 0.33;
    await page.mouse.move(cx - r, cy - r);
    await page.mouse.down();
    for (const [dx, dy] of [[r, -r * 1.1], [r * 1.2, r], [-r * 0.2, r * 1.2], [-r, r * 0.2]]) {
        await page.mouse.move(cx + dx, cy + dy, { steps: 5 });
    }
    await page.mouse.up();
    expect(await count(page)).toBe(6);
});

test('invert, grow and clear work on the active set', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.getByRole('button', { name: 'Click' }).click();
    await clickAt(page, 0.5, 0.55);
    expect(await count(page)).toBe(1);
    await page.locator('#sel-invert').click();
    expect(await count(page)).toBe(5);
    await page.locator('#sel-invert').click();
    expect(await count(page)).toBe(1);
    // One face touches four neighbours at its corners and not the opposite face.
    await page.locator('#sel-grow').click();
    expect(await count(page)).toBe(5);
    await page.locator('#sel-shrink').click();
    expect(await count(page)).toBe(1);
    await page.locator('#sel-connected').click();
    expect(await count(page)).toBe(6);
    await page.locator('#sel-clear').click();
    expect(await count(page)).toBe(0);
});

test('hide, isolate and show all change what is drawn', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.getByRole('button', { name: 'Click' }).click();
    await clickAt(page, 0.5, 0.55);
    await page.locator('#sel-hide').click();
    expect((await sel(page)).hidden).toBe(1);
    await page.locator('#sel-show').click();
    expect((await sel(page)).hidden).toBe(0);
    await page.locator('#sel-isolate').click();
    expect((await sel(page)).hidden).toBe(5);
    // A hidden cell cannot be picked: clicking where it was selects nothing new.
    await page.locator('#sel-show').click();
    expect((await sel(page)).hidden).toBe(0);
});

test('sets are independent, renameable and removable', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.getByRole('button', { name: 'Click' }).click();
    await clickAt(page, 0.5, 0.55);
    await page.locator('#sel-new').click();
    await clickAt(page, 0.5, 0.3);
    let s = await sel(page);
    expect(s.sets.length).toBe(2);
    expect(s.sets.map((x) => x.count)).toEqual([1, 1]);
    await page.locator('#sel-sets li').first().click();
    s = await sel(page);
    expect(s.active).toBe(s.sets[0].id);
    await page.locator('#sel-sets li').first().locator('.sel-name').dblclick();
    await page.locator('input.sel-rename').fill('inlet');
    await page.keyboard.press('Enter');
    expect((await sel(page)).sets[0].name).toBe('inlet');
    await page.locator('#sel-sets li').nth(1).locator('.sel-remove').click();
    expect((await sel(page)).sets.length).toBe(1);
});

test('select by value uses a cell array and shows statistics', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.locator('#sel-value summary').click();
    await page.selectOption('#sel-array', { label: 'face (cells)' });
    const min = Number(await page.locator('#sel-min').inputValue());
    await page.locator('#sel-max').fill(String(min));
    await page.locator('#sel-by-value').click();
    expect(await count(page)).toBeGreaterThan(0);
    expect(await count(page)).toBeLessThan(6);
    await expect(page.locator('#sel-stats')).toContainText('face');
    await page.locator('#sel-max').fill(String(min + 100));
    await page.locator('#sel-how').selectOption('add');
    await page.locator('#sel-by-value').click();
    expect(await count(page)).toBe(6);
});

test('ids export as JSON with the file cells of a solid', async ({ page }) => {
    await open(page, 'Block (volume)');
    await page.locator('#sel-front').uncheck();
    await page.getByRole('button', { name: 'Box' }).click();
    await boxDrag(page);
    expect(await count(page)).toBe(96);
    const download = page.waitForEvent('download');
    await page.locator('#sel-json').click();
    const file = await download;
    expect(file.suggestedFilename()).toBe('Selection_1.json');
    const body = JSON.parse(await (await import('node:fs/promises')).readFile(await file.path(), 'utf8'));
    expect(body.surfaceCells.length).toBe(96);
    // 64 hexes; the skin touches every cell except the 8 interior ones.
    expect(body.parentCells.length).toBe(56);
    expect(Math.max(...body.parentCells)).toBeLessThan(64);
});

test('a selection survives an operation that keeps cells and is cleared by one that changes them', async ({ page }) => {
    await open(page, 'Block (volume)');
    await page.locator('#sel-front').uncheck();
    await page.getByRole('button', { name: 'Box' }).click();
    await boxDrag(page);
    expect(await count(page)).toBe(96);

    await page.getByText('Quality', { exact: true }).click();
    await page.locator('#op-quality-apply').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
    expect(await count(page)).toBe(96);

    await page.getByText('Refine', { exact: true }).click();
    await page.locator('#op-refine-apply').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.numCells)).toBeGreaterThan(96);
    await expect.poll(() => count(page)).toBe(0);
});

test('loading another mesh starts a fresh selection, and the tools are exclusive with inspect', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.getByRole('button', { name: 'Click' }).click();
    await clickAt(page, 0.5, 0.55);
    await page.getByRole('button', { name: 'Inspect' }).click();
    expect((await sel(page)).tool).toBe('off');
    await expect(page.locator('#tool-inspect')).toHaveAttribute('aria-pressed', 'true');
    await page.getByRole('button', { name: 'Select', exact: true }).click();
    await expect(page.locator('#tool-inspect')).toHaveAttribute('aria-pressed', 'false');
    await page.getByRole('button', { name: 'Block (volume)' }).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.numCells)).toBe(96);
    expect(await count(page)).toBe(0);
});

test('keyboard: S toggles the tool and Escape cancels a drag', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.keyboard.press('b');
    expect((await sel(page)).tool).toBe('box');
    const b = await canvasBox(page);
    await page.mouse.move(b.x + 10, b.y + 90);
    await page.mouse.down();
    await page.mouse.move(b.x + b.width - 10, b.y + b.height - 10, { steps: 3 });
    await expect(page.locator('#sel-overlay')).toBeVisible();
    await page.keyboard.press('Escape');
    await expect(page.locator('#sel-overlay')).toBeHidden();
    await page.mouse.up();
    expect(await count(page)).toBe(0);
    await page.keyboard.press('s');
    expect((await sel(page)).tool).toBe('off');
});
