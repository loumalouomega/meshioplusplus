/**
 * The companion-process client (`meshioplusplus-mcp --http`, doc/dashboard.md).
 *
 * The server is optional: the page works with no server at all (the WASM
 * path), and every server-dependent control stays hidden until a `health`
 * probe succeeds. What the server adds is a second, full-strength producer
 * of the same summaries the browser computes for itself, and — later — the
 * training features that cannot run client-side.
 *
 * `fetch` is injectable so the client is unit-tested without a network
 * (`tests/viewer/unit/api.test.mjs`); the page passes nothing.
 */

export interface ServerHealth {
    version: string;
    root: string | null;
    runs_dir: string | null;
    tools: string[];
    mcp: string;
    transport: string;
    auth: 'token' | 'none';
}

/** One manifest the server found (`dataset_find`). */
export interface ServerManifest {
    path: string;
    relpath: string;
    sha256: string;
    name: string | null;
    num_entries: number;
    splits: Record<string, number>;
    mtime: number;
}

/** A tool failure — a payload, never an HTTP error (the MCP rule). */
export interface ToolError {
    error: string;
    error_type: string;
}

export function isToolError(value: unknown): value is ToolError {
    return (
        typeof value === 'object' &&
        value !== null &&
        typeof (value as ToolError).error === 'string' &&
        typeof (value as ToolError).error_type === 'string'
    );
}

export type FetchLike = (input: string, init?: RequestInit) => Promise<Response>;

const HEALTH_TIMEOUT_MS = 5_000;
/** A health scan stages every entry of a manifest; give it room. */
const TOOL_TIMEOUT_MS = 600_000;

export class ApiClient {
    readonly baseUrl: string;

    constructor(
        baseUrl: string,
        readonly token: string | null,
        private readonly fetchImpl: FetchLike = (input, init) => fetch(input, init),
    ) {
        this.baseUrl = baseUrl.replace(/\/+$/, '');
    }

    private headers(json: boolean): Record<string, string> {
        const headers: Record<string, string> = {};
        if (this.token) headers.Authorization = `Bearer ${this.token}`;
        if (json) headers['Content-Type'] = 'application/json';
        return headers;
    }

    private async request(path: string, init: RequestInit, timeoutMs: number): Promise<unknown> {
        const controller = new AbortController();
        const timer = setTimeout(() => controller.abort(), timeoutMs);
        let response: Response;
        try {
            response = await this.fetchImpl(`${this.baseUrl}${path}`, {
                ...init,
                signal: controller.signal,
            });
        } catch (e) {
            throw new Error(
                `meshio++: cannot reach the companion process at ${this.baseUrl}: ` +
                    (e instanceof Error ? e.message : String(e)),
            );
        } finally {
            clearTimeout(timer);
        }
        let body: unknown = null;
        try {
            body = await response.json();
        } catch {
            body = null;
        }
        if (!response.ok) {
            const message = isToolError(body) ? body.error : `HTTP ${response.status}`;
            throw new Error(`meshio++: companion process: ${message}`);
        }
        return body;
    }

    /** The connect probe. */
    async health(): Promise<ServerHealth> {
        const body = await this.request('/api/health', { method: 'GET', headers: this.headers(false) }, HEALTH_TIMEOUT_MS);
        if (!body || typeof body !== 'object' || !Array.isArray((body as ServerHealth).tools)) {
            throw new Error('meshio++: the companion process answered, but not with a health report');
        }
        return body as ServerHealth;
    }

    /** Call a registry tool; a tool failure comes back as a `ToolError`. */
    async tool<T = Record<string, unknown>>(
        name: string,
        args: Record<string, unknown> = {},
    ): Promise<T | ToolError> {
        const body = await this.request(
            `/api/tools/${encodeURIComponent(name)}`,
            { method: 'POST', headers: this.headers(true), body: JSON.stringify(args) },
            TOOL_TIMEOUT_MS,
        );
        return body as T | ToolError;
    }

    /** A download URL for a sandboxed file (the one place the token rides
     * the query string — an `<a download>` cannot set a header). */
    fileUrl(path: string): string {
        const query = new URLSearchParams({ path });
        if (this.token) query.set('token', this.token);
        return `${this.baseUrl}/api/files?${query.toString()}`;
    }
}
