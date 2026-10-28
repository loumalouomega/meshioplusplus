import assert from 'node:assert/strict';
import { test } from 'node:test';

import { VIEW_KEYS, VIEW_NAMES, isCameraState, parseBookmarks, parseCameraJson, viewPose, withBookmark } from '../../../src/viewer/src/render/camera.ts';
import { bandCount, bandPosition, canLogScale, legendTicks } from '../../../src/viewer/src/render/scale.ts';

const good = { position: [1, 2, 3], focalPoint: [0, 0, 0], viewUp: [0, 0, 1], parallel: false };

test('every view has a pose with a usable up vector, and the digit keys cover them', () => {
    for (const name of VIEW_NAMES) {
        const { direction, up } = viewPose(name);
        assert.ok(Math.hypot(...direction) > 0);
        // Looking straight along an axis must not leave "up" parallel to the view.
        const d = Math.hypot(...direction);
        const dot = Math.abs(direction[0] * up[0] + direction[1] * up[1] + direction[2] * up[2]) / d / Math.hypot(...up);
        if (name !== 'iso') assert.ok(dot < 1e-9, `${name}: up is perpendicular to the view`);
        else assert.ok(dot < 0.99);
    }
    assert.deepEqual(Object.values(VIEW_KEYS).sort(), [...VIEW_NAMES].sort());
});

test('standard views follow the Z-up convention', () => {
    assert.deepEqual(viewPose('front').direction, [0, -1, 0]);
    assert.deepEqual(viewPose('top').direction, [0, 0, 1]);
    assert.deepEqual(viewPose('right').direction, [1, 0, 0]);
});

test('camera states are validated before they touch the camera', () => {
    assert.ok(isCameraState(good));
    assert.ok(isCameraState({ ...good, parallel: true, parallelScale: 2 }));
    assert.ok(!isCameraState(null));
    assert.ok(!isCameraState({ ...good, position: [0, 0, 0] }), 'eye on the focal point');
    assert.ok(!isCameraState({ ...good, viewUp: [0, 0, 0] }));
    assert.ok(!isCameraState({ ...good, position: [1, 2, NaN] }));
    assert.ok(!isCameraState({ ...good, parallel: 'yes' }));
    assert.ok(!isCameraState({ ...good, parallelScale: -1 }));
    assert.ok(!isCameraState({ ...good, viewAngle: 0 }));
});

test('pasted JSON that is not a camera is rejected, not guessed at', () => {
    assert.deepEqual(parseCameraJson(JSON.stringify(good)), good);
    assert.equal(parseCameraJson('not json'), null);
    assert.equal(parseCameraJson('{"position":[1,2,3]}'), null);
});

test('bookmarks keep the valid entries and replace by name', () => {
    const list = parseBookmarks([{ name: 'a', camera: good }, { name: ' ', camera: good }, { name: 'b' }, 7, null]);
    assert.deepEqual(list.map((b) => b.name), ['a']);
    const next = withBookmark(list, { name: 'a', camera: { ...good, parallel: true } });
    assert.equal(next.length, 1);
    assert.equal(next[0].camera.parallel, true);
    assert.equal(withBookmark(next, { name: 'z', camera: good }).length, 2);
    assert.deepEqual(parseBookmarks('nope'), []);
});

test('log scaling needs a positive range and spaces ticks by decades', () => {
    assert.ok(canLogScale(1, 1000));
    assert.ok(!canLogScale(0, 10));
    assert.ok(!canLogScale(-1, 10));
    assert.ok(!canLogScale(5, 5));
    const t = legendTicks(1, 1000, 4, true);
    assert.ok(Math.abs(t[1] - 10) < 1e-9 && Math.abs(t[2] - 100) < 1e-9);
    assert.deepEqual(legendTicks(0, 10, 3, false), [0, 5, 10]);
    // A log request on an impossible range falls back to linear rather than NaN.
    assert.deepEqual(legendTicks(-1, 1, 3, true), [-1, 0, 1]);
});

test('bands: counts and positions', () => {
    assert.equal(bandCount('0'), 0);
    assert.equal(bandCount('5'), 5);
    assert.equal(bandCount('1'), 0);
    assert.equal(bandCount('abc'), 0);
    assert.equal(bandPosition(0.3, 0), 0.3);
    assert.equal(bandPosition(0, 4), 0.125);
    assert.equal(bandPosition(1, 4), 0.875);
    assert.equal(bandPosition(0.5, 2), 0.75);
});
