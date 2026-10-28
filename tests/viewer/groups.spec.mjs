/**
 * Groups: the file's named regions and cell blocks mapped onto the drawn
 * surface, plus groups saved from a selection.
 *
 * `zones.msh` is a 4x2x2 hexahedral box (16 hexes, 40 boundary quads) with two
 * volume zones, and an inlet and an outlet patch that are separate quad cells
 * in the file -- so the inlet/outlet groups only reach the surface through the
 * node-key match, not through an owner cell.
 */
import { expect, test } from '@playwright/test';

test.use({ serviceWorkers: 'block' });

async function open(page, label) {
    await page.goto('./');
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('idle');
    await page.getByRole('button', { name: label }).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
}

const groups = (page) => page.evaluate(() => window.__viewerState.groups);
const sel = (page) => page.evaluate(() => window.__viewerState.selection);
const byName = (state, name) => state.groups.find((g) => g.name === name);

async function openZones(page) {
    await open(page, 'Zones (regions)');
    await expect.poll(async () => (await groups(page)).status).toBe('ready');
}

test('regions and blocks are mapped onto the skin, patches through their nodes', async ({ page }) => {
    await openZones(page);
    const g = await groups(page);
    expect(await page.evaluate(() => window.__viewerState.numCells)).toBe(40);
    expect(byName(g, 'zone_left').cells).toBe(20);
    expect(byName(g, 'zone_right').cells).toBe(20);
    expect(byName(g, 'inlet').cells).toBe(4);
    expect(byName(g, 'outlet').cells).toBe(4);
    const blocks = g.groups.filter((x) => x.list === 'blocks');
    expect(blocks.map((x) => x.name).sort()).toEqual(['hexahedron', 'quad']);
    expect(byName(g, 'hexahedron').cells).toBe(40);
    expect(byName(g, 'quad').cells).toBe(8);
    await expect(page.locator('#groups-regions li')).toHaveCount(4);
    await expect(page.locator('#groups-blocks li')).toHaveCount(2);
});

test('hiding a group hides its cells; isolate shows only that group; show all restores', async ({ page }) => {
    await openZones(page);
    await page.getByLabel('Show zone_left').uncheck();
    expect((await sel(page)).hidden).toBe(20);
    expect(byName(await groups(page), 'zone_left').visible).toBe(false);
    await page.getByLabel('Show zone_left').check();
    expect((await sel(page)).hidden).toBe(0);

    await page.locator('#groups-regions li', { hasText: 'zone_right' }).hover();
    await page.locator('#groups-regions li', { hasText: 'zone_right' }).getByTitle('Show only this group').click();
    expect((await sel(page)).hidden).toBe(20);
    await page.locator('#groups-show-all').click();
    expect((await sel(page)).hidden).toBe(0);
});

test('a group selects its cells and frames them', async ({ page }) => {
    await openZones(page);
    const row = page.locator('#groups-regions li', { hasText: 'inlet' });
    await row.hover();
    await row.getByTitle('Select its cells').click();
    const s = await sel(page);
    expect(s.sets.find((x) => x.id === s.active).count).toBe(4);
    await row.getByTitle('Zoom to this group').click();
});

test('colour by group shows a categorical legend with the group names', async ({ page }) => {
    await openZones(page);
    await page.selectOption('#color-by', { label: 'Groups: regions' });
    await expect(page.locator('#legend')).toBeVisible();
    await expect(page.locator('#legend')).toHaveClass(/categorical/);
    await expect(page.locator('#legend .legend-cats')).toContainText('zone_left');
    await expect(page.locator('#legend .legend-cats')).toContainText('inlet');
    // Going back to a plain colouring restores the gradient legend rules.
    await page.selectOption('#color-by', { label: 'Solid colour' });
    await expect(page.locator('#legend')).toBeHidden();
});

test('a selection can be kept as a group of its own', async ({ page }) => {
    await openZones(page);
    await page.locator('#sel-front').uncheck();
    await page.getByRole('button', { name: 'Box' }).click();
    const b = await page.locator('#render').boundingBox();
    await page.mouse.move(b.x + 10, b.y + 90);
    await page.mouse.down();
    await page.mouse.move(b.x + b.width - 10, b.y + b.height - 10, { steps: 4 });
    await page.mouse.up();
    await page.locator('#sel-group').click();
    const g = await groups(page);
    const mine = g.groups.find((x) => x.user);
    expect(mine.name).toBe('Selection 1');
    expect(mine.cells).toBe(40);
    await expect(page.locator('#groups-regions li')).toHaveCount(5);
});

test('groups follow operations that keep cells and step aside for those that change them', async ({ page }) => {
    await openZones(page);
    await page.getByText('Quality', { exact: true }).click();
    await page.locator('#op-quality-apply').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
    expect((await groups(page)).status).toBe('ready');
    expect(byName(await groups(page), 'zone_left').cells).toBe(20);

    await page.getByText('Refine', { exact: true }).click();
    await page.locator('#op-refine-apply').click();
    await expect.poll(async () => (await groups(page)).status).toBe('unavailable');
    await expect(page.locator('#groups-note')).toContainText('as opened');
    expect((await groups(page)).groups.filter((x) => !x.user).length).toBe(0);

    // Undo the refine: the quality step is still there, the groups come back.
    await page.locator('#ops-undo').click();
    await expect.poll(async () => (await groups(page)).status).toBe('ready');
    expect(byName(await groups(page), 'zone_right').cells).toBe(20);
});

test('a solid with no regions still lists its cell blocks', async ({ page }) => {
    await open(page, 'Block (volume)');
    await expect.poll(async () => (await groups(page)).status).toBe('ready');
    const g = await groups(page);
    expect(g.groups.filter((x) => x.list === 'regions').length).toBe(0);
    expect(byName(g, 'hexahedron').cells).toBe(96);
    await expect(page.locator('#groups-note')).toContainText('no named regions');
});

test('a surface mesh maps its blocks directly', async ({ page }) => {
    await open(page, 'Cube (surface)');
    await expect.poll(async () => (await groups(page)).status).toBe('ready');
    const g = await groups(page);
    expect(byName(g, 'quad').cells).toBe(6);
});

test('loading another mesh replaces the groups', async ({ page }) => {
    await openZones(page);
    await page.getByRole('button', { name: 'Cube (surface)' }).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.numCells)).toBe(6);
    await expect.poll(async () => byName(await groups(page), 'zone_left')).toBeUndefined();
});

test('refine can be limited to a file region or to the selection', async ({ page }) => {
    await openZones(page);
    await page.getByText('Refine', { exact: true }).click();
    await page.locator('#op-refine-target').focus();
    await expect(page.locator('#op-refine-target option', { hasText: 'Region zone_left' })).toHaveCount(1);
    await page.selectOption('#op-refine-target', { label: 'Region zone_left' });
    await expect(page.locator('#op-refine-array')).toBeDisabled();
    await page.locator('#op-refine-apply').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.numCells)).toBeGreaterThan(40);
    // Half the mesh refined, not all of it: fewer cells than the uniform 8-for-1.
    expect(await page.evaluate(() => window.__viewerState.numCells)).toBeLessThan(40 * 4);
    await expect(page.locator('.op-chip').first()).toContainText('refine · region zone_left');
});

test('refining the selection sends its file cells', async ({ page }) => {
    await openZones(page);
    const row = page.locator('#groups-regions li', { hasText: 'zone_right' });
    await row.hover();
    await row.getByTitle('Select its cells').click();
    await page.getByText('Refine', { exact: true }).click();
    await page.locator('#op-refine-target').focus();
    await expect(page.locator('#op-refine-target option', { hasText: /^Selection \(\d+ cells\)$/ })).toHaveCount(1);
    await page.selectOption('#op-refine-target', { value: 'selection' });
    await page.locator('#op-refine-apply').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.numCells)).toBeGreaterThan(40);
    await expect(page.locator('.op-chip').first()).toContainText('refine · 8 cells');
});

test('saved groups are written with the mesh, and the result says what was kept', async ({ page }) => {
    await openZones(page);
    await page.locator('#convert-groups-row').waitFor({ state: 'hidden' });
    const row = page.locator('#groups-regions li', { hasText: 'zone_left' });
    await row.hover();
    await row.getByTitle('Select its cells').click();
    await page.locator('#sel-group').click();
    await expect(page.locator('#convert-groups-row')).toBeVisible();

    await page.selectOption('#convert-format', 'gmsh');
    const download = page.waitForEvent('download');
    await page.getByRole('button', { name: 'Download' }).click();
    const file = await download;
    const text = await (await import('node:fs/promises')).readFile(await file.path(), 'latin1');
    expect(file.suggestedFilename()).toBe('zones.gmsh');
    // gmsh 4.1 is binary but names are ASCII; the saved group is among them.
    expect(text).toContain('Selection 1');
    await expect(page.locator('#status')).toContainText('wrote 1 group');
});

test('the saved-groups option disappears once an operation is applied', async ({ page }) => {
    await openZones(page);
    await page.locator('#sel-front').uncheck();
    await page.getByRole('button', { name: 'Box' }).click();
    const b = await page.locator('#render').boundingBox();
    await page.mouse.move(b.x + 10, b.y + 90);
    await page.mouse.down();
    await page.mouse.move(b.x + b.width - 10, b.y + b.height - 10, { steps: 3 });
    await page.mouse.up();
    await page.locator('#sel-group').click();
    await expect(page.locator('#convert-groups-row')).toBeVisible();
    await page.getByText('Quality', { exact: true }).click();
    await page.locator('#op-quality-apply').click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
    await expect(page.locator('#convert-groups-row')).toBeHidden();
});
