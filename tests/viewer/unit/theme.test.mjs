import assert from 'node:assert/strict';
import { test } from 'node:test';

import { SCENE_THEMES, THEME_MODES, cssColor, isThemeMode, resolveTheme } from '../../../src/viewer/src/ui/theme.ts';

test('auto follows the OS preference, explicit modes ignore it', () => {
    assert.equal(resolveTheme('auto', true), 'dark');
    assert.equal(resolveTheme('auto', false), 'light');
    assert.equal(resolveTheme('dark', false), 'dark');
    assert.equal(resolveTheme('light', true), 'light');
    assert.equal(resolveTheme('scientific', true), 'scientific');
});

test('isThemeMode accepts only the four modes', () => {
    for (const mode of THEME_MODES) assert.ok(isThemeMode(mode));
    assert.ok(!isThemeMode('solarized'));
    assert.ok(!isThemeMode(null));
});

test('every scene theme has a background, a body colour and eight palette entries', () => {
    for (const theme of Object.values(SCENE_THEMES)) {
        assert.equal(theme.background.length, 3);
        assert.equal(theme.body.length, 3);
        assert.equal(theme.palette.length, 8);
        for (const c of [...theme.background, ...theme.body, ...theme.palette.flat()]) {
            assert.ok(c >= 0 && c <= 1);
        }
    }
});

test('the scientific scene is white and the light scene is the extension grey', () => {
    assert.deepEqual(SCENE_THEMES.scientific.background, [1, 1, 1]);
    assert.equal(cssColor(SCENE_THEMES.light.background), 'rgb(240, 240, 240)');
    assert.equal(cssColor(SCENE_THEMES.dark.background), 'rgb(30, 30, 30)');
});

test('cssColor clamps out-of-range channels', () => {
    assert.equal(cssColor([2, -1, 0.5]), 'rgb(255, 0, 128)');
});
