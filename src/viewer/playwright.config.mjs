import { defineConfig, devices } from '@playwright/test';

// The base path is load-bearing: `vite preview` serves the app under the same
// prefix it will have on GitHub Pages, and a bare `/` 404s.
const BASE = '/meshioplusplus/viewer/';

export default defineConfig({
    testDir: '../../tests/viewer',
    // Only the browser specs. Playwright's default pattern also matches
    // `*.test.mjs`, which would drag in the `node --test` unit tests under
    // tests/viewer/unit/ — they import TypeScript directly and are run by
    // `npm run test:unit` instead.
    testMatch: '**/*.spec.mjs',
    // The WASM module has to instantiate before the first mesh can be read,
    // which is slower than a typical DOM assertion on a cold CI runner.
    timeout: 60_000,
    expect: { timeout: 20_000 },
    forbidOnly: !!process.env.CI,
    retries: process.env.CI ? 1 : 0,
    reporter: process.env.CI ? [['list'], ['html', { open: 'never' }]] : 'list',
    use: {
        baseURL: `http://localhost:4173${BASE}`,
        trace: 'retain-on-failure',
    },
    projects: [{ name: 'chromium', use: { ...devices['Desktop Chrome'] } }],
    webServer: {
        command: 'npm run preview',
        url: `http://localhost:4173${BASE}`,
        reuseExistingServer: !process.env.CI,
        timeout: 60_000,
    },
});
