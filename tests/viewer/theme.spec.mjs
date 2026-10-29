/**
 * Look and feel: both colour schemes, against the design tokens.
 *
 * The viewer shares its metric scale and its VS Code dark/light colours with the
 * Kratos MDPA Preview extension (src/viewer/src/tokens.css). These assertions
 * read computed styles, so they need no GPU.
 */
import { expect, test } from '@playwright/test';

test.use({ serviceWorkers: 'block' });

const style = (page, selector, prop) =>
    page.locator(selector).first().evaluate((el, p) => getComputedStyle(el)[p], prop);

async function idle(page) {
    await page.goto('./');
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('idle');
}

for (const [scheme, panel, header] of [
    ['dark', 'rgb(37, 37, 38)', 'rgb(37, 37, 38)'],
    ['light', 'rgb(243, 243, 243)', 'rgb(243, 243, 243)'],
]) {
    test(`auto follows the OS: ${scheme}`, async ({ page }) => {
        await page.emulateMedia({ colorScheme: scheme });
        await idle(page);
        await expect(page.locator('html')).toHaveAttribute('data-theme', scheme);
        expect(await style(page, '#panel', 'backgroundColor')).toBe(panel);
        expect(await style(page, '#menubar', 'backgroundColor')).toBe(header);
    });
}

test('tokens: radii and type sizes match the extension scale', async ({ page }) => {
    await page.emulateMedia({ colorScheme: 'dark' });
    await idle(page);
    // --ui-r-sm for buttons and controls, --ui-r-md for floating cards.
    expect(await style(page, '.sample', 'borderRadius')).toBe('5px');
    expect(await style(page, '#viewport-tools', 'borderRadius')).toBe('8px');
    // --ui-fs-xl body text, --ui-fs-md for hints.
    expect(await style(page, 'body', 'fontSize')).toBe('14px');
    expect(await style(page, '.hint', 'fontSize')).toBe('12px');
    expect(await style(page, '#statusbar', 'fontSize')).toBe('11px');
});

test('the scheme picker overrides the OS and survives a reload', async ({ page }) => {
    await page.emulateMedia({ colorScheme: 'dark' });
    await idle(page);
    await page.selectOption('#theme-select', 'light');
    await expect(page.locator('html')).toHaveAttribute('data-theme', 'light');
    expect(await style(page, '#panel', 'backgroundColor')).toBe('rgb(243, 243, 243)');
    await page.reload();
    await expect(page.locator('html')).toHaveAttribute('data-theme', 'light');
    await expect(page.locator('#theme-select')).toHaveValue('light');
    await page.selectOption('#theme-select', 'scientific');
    await expect(page.locator('html')).toHaveAttribute('data-theme', 'scientific');
});

test('menus open, check items reflect state and Escape closes', async ({ page }) => {
    await idle(page);
    await page.getByRole('button', { name: 'View' }).click();
    await expect(page.locator('#menu-view')).toBeVisible();
    const edges = page.getByRole('menuitemcheckbox', { name: 'Show edges' });
    await expect(edges).toHaveAttribute('aria-checked', 'false');
    await edges.click();
    await expect(page.locator('#edges')).toBeChecked();
    // A check item keeps the menu open, with its tick updated.
    await expect(page.getByRole('menuitemcheckbox', { name: 'Show edges' })).toHaveAttribute('aria-checked', 'true');
    await page.keyboard.press('Escape');
    await expect(page.locator('#menu-view')).toBeHidden();
});

test('sections collapse with aria-expanded and remember it', async ({ page }) => {
    await idle(page);
    await page.getByRole('button', { name: 'Cube (surface)' }).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
    const chevron = page.getByRole('button', { name: 'Toggle Operations' });
    await expect(chevron).toHaveAttribute('aria-expanded', 'true');
    await chevron.click();
    await expect(chevron).toHaveAttribute('aria-expanded', 'false');
    await expect(page.locator('#ops-section .sb-body')).toBeHidden();
    await page.reload();
    await page.getByRole('button', { name: 'Cube (surface)' }).click();
    await expect(page.getByRole('button', { name: 'Toggle Operations' })).toHaveAttribute('aria-expanded', 'false');
});

test('keyboard shortcuts act and the sheet lists them', async ({ page }) => {
    await idle(page);
    await page.keyboard.press('?');
    await expect(page.locator('#shortcuts-dialog')).toBeVisible();
    await expect(page.locator('#shortcuts-body')).toContainText('Fit the mesh to the window');
    await page.keyboard.press('Escape');
    await expect(page.locator('#shortcuts-dialog')).toBeHidden();
    await page.keyboard.press('e');
    await expect(page.locator('#edges')).toBeChecked();
});

for (const scheme of ['dark', 'light']) {
    test(`the dataset manager loads cleanly in ${scheme}`, async ({ page }) => {
        const errors = [];
        page.on('pageerror', (e) => errors.push(e.message));
        await page.emulateMedia({ colorScheme: scheme });
        await page.goto('./dataset.html');
        await expect(page.locator('#ws-pick')).toBeVisible();
        // The manager never sets data-theme, so it stays on the dark palette.
        expect(await style(page, '#panel', 'backgroundColor')).toBe('rgb(37, 37, 38)');
        expect(errors.filter((m) => !/scope/.test(m))).toEqual([]);
    });
}

test('the status bar reports what is open', async ({ page }) => {
    await idle(page);
    await page.getByRole('button', { name: 'Cube (surface)' }).click();
    await expect.poll(() => page.evaluate(() => window.__viewerState.status)).toBe('ready');
    await expect(page.locator('#sb-file')).toHaveText('cube.vtp');
    await expect(page.locator('#sb-points')).toContainText('8 points');
    await expect(page.locator('#sb-cells')).toContainText('6 surface cells');
});
