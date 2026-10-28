import assert from 'node:assert/strict';
import { test } from 'node:test';

import { Store } from '../../../src/viewer/src/app/store.ts';
import { findByShortcut, matchShortcut, registerAction, runAction } from '../../../src/viewer/src/app/actions.ts';
import { clampSidebarWidth } from '../../../src/viewer/src/features/shell.ts';

const key = (k, mods = {}) => ({ key: k, ctrlKey: false, metaKey: false, altKey: false, shiftKey: false, ...mods });

test('a store notifies only on change and can be unsubscribed', () => {
    const store = new Store(1);
    const seen = [];
    const off = store.subscribe((v, p) => seen.push([v, p]));
    store.set(1);
    store.set(2);
    store.update((v) => v + 1);
    off();
    store.set(9);
    assert.deepEqual(seen, [[2, 1], [3, 2]]);
    assert.equal(store.get(), 9);
});

test('shortcuts match letters, digits, named keys and modifiers', () => {
    assert.ok(matchShortcut('F', key('f')));
    assert.ok(matchShortcut('1', key('1')));
    assert.ok(matchShortcut('Esc', key('Escape')));
    assert.ok(matchShortcut('Ctrl+O', key('o', { ctrlKey: true })));
    assert.ok(matchShortcut('Ctrl+O', key('o', { metaKey: true })));
    assert.ok(!matchShortcut('Ctrl+O', key('o')));
    assert.ok(!matchShortcut('F', key('f', { ctrlKey: true })));
    assert.ok(matchShortcut('Shift+I', key('I', { shiftKey: true })));
    assert.ok(matchShortcut('?', key('?', { shiftKey: true })));
});

test('actions run through the registry and honour enabled()', () => {
    let ran = 0;
    let on = false;
    registerAction({ id: 't.a', label: 'A', shortcut: 'Q', enabled: () => on, run: () => ran++ });
    assert.equal(runAction('t.a'), false);
    assert.equal(findByShortcut(key('q')), undefined);
    on = true;
    assert.equal(findByShortcut(key('q'))?.id, 't.a');
    assert.equal(runAction('t.a'), true);
    assert.equal(ran, 1);
    assert.equal(runAction('t.missing'), false);
});

test('the sidebar width is clamped and leaves room for the viewport', () => {
    assert.equal(clampSidebarWidth(50), 220);
    assert.equal(clampSidebarWidth(2000), 640);
    assert.equal(clampSidebarWidth(500, 600), 360);
    assert.equal(clampSidebarWidth(300.4), 300);
});
