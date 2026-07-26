/*! coi-serviceworker v0.1.7 - Guido Zuidhof and contributors, licensed under MIT */
// Makes a statically-hosted page cross-origin isolated by installing a service
// worker that re-serves every response with the COOP/COEP headers
// SharedArrayBuffer requires. This is what lets the meshio++ browser viewer load
// the THREADED wasm artifact (meshioplusplus_wasm_mt) on GitHub Pages, which
// cannot itself set those headers. Without this the page is not isolated and the
// loader (src/wasm/src/index.mjs) auto-falls-back to the sequential build, so
// the SW is a performance enhancement, never a correctness requirement.
//
// It is deliberately vendored verbatim (upstream: gzuidhof/coi-serviceworker)
// rather than added as a dependency: a single-file public-domain-ish MIT snippet
// is easier to audit than a build-time package, and it is served only by the
// web build (src/viewer's `public/` is excluded from the single-file embed
// build, which never touches wasm at all).
let coepCredentialless = false;
if (typeof window === 'undefined') {
    self.addEventListener('install', () => self.skipWaiting());
    self.addEventListener('activate', (event) => event.waitUntil(self.clients.claim()));

    self.addEventListener('message', (ev) => {
        if (!ev.data) {
            return;
        } else if (ev.data.type === 'deregister') {
            self.registration
                .unregister()
                .then(() => self.clients.matchAll())
                .then((clients) => {
                    clients.forEach((client) => client.navigate(client.url));
                });
        } else if (ev.data.type === 'coepCredentialless') {
            coepCredentialless = ev.data.value;
        }
    });

    self.addEventListener('fetch', function (event) {
        const r = event.request;
        if (r.cache === 'only-if-cached' && r.mode !== 'same-origin') {
            return;
        }

        const request =
            coepCredentialless && r.mode === 'no-cors'
                ? new Request(r, {
                      credentials: 'omit',
                  })
                : r;
        event.respondWith(
            fetch(request)
                .then((response) => {
                    if (response.status === 0) {
                        return response;
                    }

                    const newHeaders = new Headers(response.headers);
                    newHeaders.set('Cross-Origin-Embedder-Policy',
                        coepCredentialless ? 'credentialless' : 'require-corp');
                    if (!coepCredentialless) {
                        newHeaders.set('Cross-Origin-Resource-Policy', 'cross-origin');
                    }
                    newHeaders.set('Cross-Origin-Opener-Policy', 'same-origin');

                    return new Response(response.body, {
                        status: response.status,
                        statusText: response.statusText,
                        headers: newHeaders,
                    });
                })
                .catch((e) => console.error(e))
        );
    });
} else {
    (() => {
        const reloadedBySelf = window.sessionStorage.getItem('coiReloadedBySelf');
        window.sessionStorage.removeItem('coiReloadedBySelf');
        const coepDegrading = reloadedBySelf == 'coepdegrade';

        // You can customize the behavior of this script through a global `coi` variable.
        const coi = {
            shouldRegister: () => !reloadedBySelf,
            shouldDeregister: () => false,
            coepCredentialless: () => true,
            coepDegrade: () => true,
            doReload: () => window.location.reload(),
            quiet: false,
            ...window.coi,
        };

        const n = navigator;
        const controlling = n.serviceWorker && n.serviceWorker.controller;

        // Record the failure if the page is served by the service worker.
        if (controlling && !window.crossOriginIsolated) {
            window.sessionStorage.setItem('coiCoepHasFailed', 'true');
        }
        const coepHasFailed = window.sessionStorage.getItem('coiCoepHasFailed');

        if (controlling) {
            // Reload only on the first failure.
            const reloadToDegrade = coi.coepDegrade() && !(coepDegrading || window.crossOriginIsolated);
            n.serviceWorker.controller.postMessage({
                type: 'coepCredentialless',
                value: (reloadToDegrade || coepHasFailed && coi.coepDegrade())
                    ? false
                    : coi.coepCredentialless(),
            });
            if (reloadToDegrade) {
                !coi.quiet && console.log('Reloading page to degrade COEP.');
                window.sessionStorage.setItem('coiReloadedBySelf', 'coepdegrade');
                coi.doReload('coepdegrade');
            }

            if (coi.shouldDeregister()) {
                n.serviceWorker.controller.postMessage({ type: 'deregister' });
            }
        }

        // If we're already coi: do nothing. Perhaps it's due to this script doing its job, or COOP/COEP are
        // already set from the origin server. Also if the browser has no notion of crossOriginIsolated, just give up here.
        if (window.crossOriginIsolated !== false || !coi.shouldRegister()) return;

        if (!window.isSecureContext) {
            !coi.quiet && console.log('COOP/COEP Service Worker not registered, a secure context is required.');
            return;
        }

        // In some environments (e.g. Firefox private mode) this won't be available
        if (!n.serviceWorker) {
            !coi.quiet &&
                console.error('COOP/COEP Service Worker not registered, perhaps due to private mode.');
            return;
        }

        n.serviceWorker.register(window.document.currentScript.src).then(
            (registration) => {
                !coi.quiet && console.log('COOP/COEP Service Worker registered', registration.scope);

                registration.addEventListener('updatefound', () => {
                    !coi.quiet && console.log('Reloading page to make use of updated COOP/COEP Service Worker.');
                    window.sessionStorage.setItem('coiReloadedBySelf', 'updatefound');
                    coi.doReload();
                });

                // If the registration is active, but it's not controlling the page
                if (registration.active && !n.serviceWorker.controller) {
                    !coi.quiet && console.log('Reloading page to make use of COOP/COEP Service Worker.');
                    window.sessionStorage.setItem('coiReloadedBySelf', 'notcontrolling');
                    coi.doReload();
                }
            },
            (err) => {
                !coi.quiet && console.error('COOP/COEP Service Worker failed to register:', err);
            }
        );
    })();
}
