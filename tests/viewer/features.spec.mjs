/**
 * The measuring mode (roadmap 7.9.1): click a feature to grow and fit it,
 * relate two fitted features, and the deviation and thickness views.
 *
 * The probe is a bore: a half-pipe trough of radius 2.5 along x with flat lips
 * at z = 0, generated here and opened as an OBJ. The top view looks down -z,
 * so the middle of the render area lands on the trough's floor. Asserts on
 * `window.__viewerState.features`, like every other spec, so no GPU is needed.
 */
import { expect, test } from '@playwright/test';

// Fitting does not need cross-origin isolation; block the COOP/COEP service
// worker so its one-off reload cannot land mid-poll.
test.use({ serviceWorkers: 'block' });

const RADIUS = 2.5;

/** The trough as OBJ text: 24 x 20 quads of half-pipe and two 4 x 20 lips. */
function troughObj() {
    const v = [];
    const f = [];
    const nx = 20;
    const grid = (rows, point) => {
        const first = v.length;
        for (let j = 0; j < rows.length; ++j)
            for (let i = 0; i <= nx; ++i) v.push(point(-5 + (10 * i) / nx, rows[j]));
        for (let j = 0; j + 1 < rows.length; ++j)
            for (let i = 0; i < nx; ++i) {
                const a = first + j * (nx + 1) + i + 1;
                f.push([a, a + 1, a + nx + 2, a + nx + 1]);
            }
    };
    const angles = Array.from({ length: 25 }, (_, k) => Math.PI + (Math.PI * k) / 24);
    grid(angles, (x, t) => [x, RADIUS * Math.cos(t), RADIUS * Math.sin(t)]);
    const lip = (from, to) => Array.from({ length: 5 }, (_, k) => from + ((to - from) * k) / 4);
    grid(lip(RADIUS, RADIUS + 2), (x, y) => [x, y, 0]);
    grid(lip(-RADIUS - 2, -RADIUS), (x, y) => [x, y, 0]);
    return [
        ...v.map((p) => `v ${p.map((c) => c.toFixed(6)).join(' ')}`),
        ...f.map((q) => `f ${q.join(' ')}`),
        '',
    ].join('\n');
}

const features = (page) => page.evaluate(() => window.__viewerState.features);
const status = (page) => page.evaluate(() => window.__viewerState.status);

async function openTrough(page) {
    await page.goto('./');
    await expect.poll(() => status(page)).toBe('idle');
    await page.locator('#file').setInputFiles({
        name: 'trough.obj',
        mimeType: 'text/plain',
        buffer: Buffer.from(troughObj()),
    });
    await expect.poll(() => status(page)).toBe('ready');
    await page.locator('#view-buttons [data-view="top"]').click();
}

/** Click at a fraction of the render area's height above its middle and wait for the fit. */
async function fitAt(page, up = 0) {
    const before = (await features(page)).elements.length;
    const box = await page.locator('#render').boundingBox();
    await page.mouse.click(box.x + box.width / 2, box.y + box.height * (0.5 - up));
    await expect
        .poll(async () => {
            const s = await features(page);
            return !s.busy && (s.elements.length > before || s.error !== null);
        })
        .toBe(true);
    return features(page);
}

test('clicking the bore fits its radius', async ({ page }) => {
    await openTrough(page);
    await expect(page.locator('#tool-fit')).toHaveAttribute('aria-pressed', 'false');
    await page.locator('#tool-fit').click();
    await expect(page.locator('#tool-fit')).toHaveAttribute('aria-pressed', 'true');
    await page.locator('#feat-type').selectOption('cylinder');

    const s = await fitAt(page);
    expect(s.error).toBeNull();
    expect(s.elements).toHaveLength(1);
    const bore = s.elements[0];
    expect(bore.kind).toBe('cylinder');
    expect(bore.radius).toBeCloseTo(RADIUS, 6);
    expect(bore.numCells).toBeGreaterThan(100);
    await expect(page.locator('#feat-list')).toContainText('Diameter');
    await expect(page.locator('#feat-list')).toContainText('5');

    // The other click tools give the click back.
    await page.locator('#tool-measure').click();
    expect((await features(page)).active).toBe(false);
});

test('a bore and a lip relate; deviation and thickness colour the surface', async ({ page }) => {
    await openTrough(page);
    await page.locator('#tool-fit').click();
    await page.locator('#feat-type').selectOption('cylinder');
    const bore = await fitAt(page);
    expect(bore.elements[0].kind).toBe('cylinder');

    // The share within the tolerance counts the trough, not the lips.
    await page.locator('#feat-views').evaluate((d) => d.setAttribute('open', ''));
    await page.locator('#feat-deviation').click();
    await expect.poll(async () => (await features(page)).view?.shareWithin ?? null).not.toBeNull();
    const dev = (await features(page)).view;
    expect(dev.kind).toBe('deviation');
    expect(dev.shareWithin).toBeGreaterThan(0.4);
    expect(dev.shareWithin).toBeLessThan(1);
    expect(await page.evaluate(() => window.__viewerState.arrays)).toContain('deviation:distance');
    await expect(page.locator('#feat-view-result')).toContainText('Within');

    await page.locator('#feat-thickness').click();
    await expect.poll(async () => (await features(page)).view?.kind).toBe('thickness');
    await expect
        .poll(() => page.evaluate(() => window.__viewerState.arrays))
        .toContain('thickness:distance');

    // Step up the screen until a click lands on a lip and fits a plane.
    await page.locator('#feat-type').selectOption('auto');
    let s = null;
    for (let k = 1; k <= 9; ++k) {
        s = await fitAt(page, 0.045 * k);
        if (s.elements.at(-1)?.kind === 'plane') break;
    }
    expect(s.elements.at(-1).kind).toBe('plane');
    // The two newest are selected; reselect the bore and the lip.
    const lip = s.elements.at(-1).name;
    for (const name of s.selected) await page.locator(`.feat-item[data-name="${name}"]`).click();
    await page.locator('.feat-item[data-name="F1"]').click();
    await page.locator(`.feat-item[data-name="${lip}"]`).click();
    const r = (await features(page)).relation;
    expect(r.what).toBe('axis to plane');
    expect(r.angleDeg).toBeCloseTo(0, 4);
    expect(r.distance).toBeCloseTo(0, 4);
    await expect(page.locator('#feat-relation')).toContainText('Angle');
});
