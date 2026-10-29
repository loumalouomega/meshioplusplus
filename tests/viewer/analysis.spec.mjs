/**
 * Data table and statistics: windowed rows, sorting, CSV, row -> selection,
 * histogram bar -> selection.
 */
import { expect, test } from '@playwright/test';

test.use({ serviceWorkers: 'block' });

async function open(page, label) {
    await page.goto('./');
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('idle');
    await page.getByRole('button', { name: label }).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
}

const analysis = (page) => page.evaluate(() => window.__viewerState.analysis);
const sel = (page) => page.evaluate(() => window.__viewerState.selection);

test('the table lists every cell, keeps the DOM small, and switches to points', async ({ page }) => {
    await open(page, 'Wave (point data)');
    await expect(page.locator('#table-panel')).toBeHidden();
    await page.click('#tool-table');
    const state = await analysis(page);
    expect(state.table.open).toBe(true);
    const cells = await page.evaluate(() => window.__viewerState.numCells);
    expect(state.table.rows).toBe(cells);
    expect(await page.locator('.tbl-row').count()).toBeLessThan(Math.min(cells, 200) + 1);

    await page.getByRole('button', { name: 'Points', exact: true }).click();
    const pts = await analysis(page);
    expect(pts.table.kind).toBe('points');
    expect(pts.table.rows).toBe(await page.evaluate(() => window.__viewerState.numPoints));
    await expect(page.locator('.tbl-th', { hasText: 'x' }).first()).toBeVisible();
});

test('clicking a header sorts, and a second click reverses', async ({ page }) => {
    await open(page, 'Wave (point data)');
    await page.click('#tool-table');
    await page.getByRole('button', { name: 'Points', exact: true }).click();
    await page.locator('.tbl-th', { hasText: /^x/ }).first().click();
    let s = (await analysis(page)).table;
    expect(s.sortedBy).toBe('x');
    expect(s.direction).toBe('asc');
    const first = await page.locator('.tbl-row').first().locator('span').nth(1).textContent();
    await page.locator('.tbl-th', { hasText: /^x/ }).first().click();
    s = (await analysis(page)).table;
    expect(s.direction).toBe('desc');
    const firstDesc = await page.locator('.tbl-row').first().locator('span').nth(1).textContent();
    expect(Number(firstDesc)).toBeGreaterThanOrEqual(Number(first));
});

test('clicking a cell row selects that cell, and CSV downloads the table', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.click('#tool-table');
    await page.locator('.tbl-row').nth(3).click();
    const s = await sel(page);
    expect(s.sets.find((x) => x.id === s.active).count).toBe(1);

    const download = page.waitForEvent('download');
    await page.click('#tbl-csv');
    const file = await download;
    expect(file.suggestedFilename()).toMatch(/-cells\.csv$/);
    const text = await (await import('node:fs/promises')).readFile(await file.path(), 'utf8');
    expect(text.split('\n')[0]).toContain('cell');
    expect(text.trim().split('\n').length).toBe(1 + (await page.evaluate(() => window.__viewerState.numCells)));
});

test('the column picker hides a column', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.click('#tool-table');
    const before = await page.locator('.tbl-th').count();
    await page.locator('#tbl-columns summary').click();
    await page.locator('#tbl-column-list input').first().uncheck();
    expect(await page.locator('.tbl-th').count()).toBe(before - 1);
});

test('the statistics histogram counts every value and a bar selects cells', async ({ page }) => {
    await open(page, 'Wave (point data)');
    await expect(page.locator('#stats-section')).toBeVisible();
    const h = (await analysis(page)).histogram;
    expect(h.array).not.toBeNull();
    expect(h.counts.length).toBe(20);
    expect(h.counts.reduce((a, b) => a + b, 0)).toBeGreaterThan(0);

    await page.selectOption('#stats-bins', '10');
    expect((await analysis(page)).histogram.counts.length).toBe(10);

    await page.locator('.hist-bar').nth(0).click();
    const s = await sel(page);
    expect(s.sets.find((x) => x.id === s.active).count).toBeGreaterThan(0);
});

test('the T key toggles the table', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await page.locator('#render').click({ position: { x: 5, y: 5 } });
    await page.keyboard.press('t');
    expect((await analysis(page)).table.open).toBe(true);
    await page.keyboard.press('t');
    expect((await analysis(page)).table.open).toBe(false);
});
