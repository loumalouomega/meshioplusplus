/**
 * The companion-process client, driven with an injected `fetch` — the
 * header it sends, the error shapes it turns into named errors, and the
 * one place the token rides a URL.
 *
 * Run with `node --test` via tests/viewer/unit/build-and-run.mjs.
 */
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { ApiClient, isToolError } from '../../../src/viewer/src/dataset/api.ts';

function fakeFetch(handler) {
    const calls = [];
    const fetchImpl = async (url, init) => {
        calls.push({ url, init });
        const { status = 200, body = null } = handler(url, init);
        return {
            ok: status >= 200 && status < 300,
            status,
            json: async () => body,
        };
    };
    return { fetchImpl, calls };
}

test('health sends the bearer token and returns the report', async () => {
    const report = { version: '10.23.0', root: '/w', runs_dir: '/w/runs', tools: ['info'], mcp: '/mcp', transport: 'streamable-http', auth: 'token' };
    const { fetchImpl, calls } = fakeFetch(() => ({ body: report }));
    const api = new ApiClient('http://127.0.0.1:8765/', 'tok', fetchImpl);
    assert.equal(api.baseUrl, 'http://127.0.0.1:8765');
    assert.deepEqual(await api.health(), report);
    assert.equal(calls[0].url, 'http://127.0.0.1:8765/api/health');
    assert.equal(calls[0].init.headers.Authorization, 'Bearer tok');
});

test('a 401 and an unreachable server become named errors', async () => {
    const denied = fakeFetch(() => ({ status: 401, body: { error: 'missing or invalid bearer token', error_type: 'PermissionError' } }));
    await assert.rejects(
        new ApiClient('http://x', 'bad', denied.fetchImpl).health(),
        /missing or invalid bearer token/,
    );
    const down = { fetchImpl: async () => { throw new TypeError('Failed to fetch'); } };
    await assert.rejects(new ApiClient('http://x', null, down.fetchImpl).health(), /cannot reach the companion process/);
    const notHealth = fakeFetch(() => ({ body: { hello: 1 } }));
    await assert.rejects(new ApiClient('http://x', null, notHealth.fetchImpl).health(), /not with a health report/);
});

test('tool calls POST JSON and pass tool errors through as payloads', async () => {
    const { fetchImpl, calls } = fakeFetch((url) =>
        url.endsWith('/api/tools/info')
            ? { body: { num_points: 5 } }
            : { body: { error: 'nope', error_type: 'ValueError' } },
    );
    const api = new ApiClient('http://x', null, fetchImpl);
    assert.deepEqual(await api.tool('info', { input_path: 'a.vtu' }), { num_points: 5 });
    assert.equal(calls[0].init.method, 'POST');
    assert.equal(calls[0].init.body, '{"input_path":"a.vtu"}');
    assert.equal(calls[0].init.headers['Content-Type'], 'application/json');
    assert.equal('Authorization' in calls[0].init.headers, false);
    const failed = await api.tool('dataset_health', {});
    assert.ok(isToolError(failed));
    assert.equal(failed.error_type, 'ValueError');
    assert.equal(isToolError({ num_points: 5 }), false);
});

test('fileUrl carries the sandboxed path and the token', () => {
    const api = new ApiClient('http://x', 't k', async () => { throw new Error('unused'); });
    assert.equal(api.fileUrl('runs/a b.mdlus'), 'http://x/api/files?path=runs%2Fa+b.mdlus&token=t+k');
    const noToken = new ApiClient('http://x', null, async () => { throw new Error('unused'); });
    assert.equal(noToken.fileUrl('a'), 'http://x/api/files?path=a');
});
