/**
 * Training-job types and the poller that follows a run.
 *
 * The companion process exposes jobs as ordinary registry tools
 * (`train_status`/`train_metrics`/`train_log`), so following a run is three
 * incremental calls on a timer: the status, the metrics rows since the last
 * epoch seen, and the log from the last byte offset. Polling stops by itself
 * once the job reaches a terminal state — a finished run must not keep a tab
 * busy.
 */

import { isToolError, type ApiClient } from './api';

export interface MetricRow {
    epoch: number;
    train_loss: number;
    valid_loss: number | null;
    lr: number;
    elapsed: number;
    epoch_seconds: number;
    timestamp: number;
}

export interface JobStatus {
    job_id: string;
    run_dir: string;
    status: string;
    pid: number | null;
    started: number | null;
    finished: number | null;
    exit_code: number | null;
    reason?: string | null;
    manifest: string | null;
    best_checkpoint: string | null;
    epoch: number | null;
    epochs: number | null;
    best_epoch: number | null;
    best_valid_loss: number | null;
    eta_seconds: number | null;
    device: string | null;
    completed: boolean;
    num_metrics: number;
    last: MetricRow | null;
}

export interface JobSummary extends JobStatus {
    fields: string[];
    target_fields: string[];
    train_split: string | null;
    valid_split: string | null;
    batch_size: number | null;
    learning_rate: number | null;
    seed: number | null;
    hidden_dim: number | null;
    processor_size: number | null;
    tags: string[];
    notes: string | null;
    final_train_loss: number | null;
    final_valid_loss: number | null;
    duration_seconds: number | null;
}

export interface CheckpointInfo {
    path: string;
    name: string;
    kind: 'periodic' | 'best' | 'final';
    epoch: number | null;
    valid_loss: number | null;
    size: number;
    is_best: boolean;
}

export const TERMINAL_STATUSES = ['finished', 'failed', 'stopped'];

export function isTerminal(status: string): boolean {
    return TERMINAL_STATUSES.includes(status);
}

export interface PollUpdate {
    status: JobStatus;
    /** Rows appended since the previous update. */
    newRows: MetricRow[];
    /** Log text appended since the previous update. */
    newLog: string;
    checkpoints: CheckpointInfo[] | null;
    terminal: boolean;
}

export interface PollerOptions {
    onUpdate(update: PollUpdate): void;
    onError?(message: string): void;
    intervalMs?: number;
}

/** Follow one job until it stops (or `stop()` is called). */
export class JobPoller {
    private timer = 0;
    private running = false;
    private sinceEpoch = 0;
    private logOffset = 0;

    constructor(
        private readonly api: ApiClient,
        readonly jobId: string,
        private readonly options: PollerOptions,
    ) {}

    start(): void {
        if (this.running) return;
        this.running = true;
        void this.tick();
    }

    stop(): void {
        this.running = false;
        if (this.timer) {
            clearTimeout(this.timer);
            this.timer = 0;
        }
    }

    private schedule(): void {
        if (!this.running) return;
        this.timer = window.setTimeout(() => void this.tick(), this.options.intervalMs ?? 2000);
    }

    private async tick(): Promise<void> {
        if (!this.running) return;
        try {
            const status = await this.api.tool<JobStatus>('train_status', {
                job_id: this.jobId,
            });
            if (isToolError(status)) throw new Error(status.error);
            const metrics = await this.api.tool<{ rows: MetricRow[] }>('train_metrics', {
                job_id: this.jobId,
                since_epoch: this.sinceEpoch,
            });
            const newRows = isToolError(metrics) ? [] : metrics.rows;
            if (newRows.length) this.sinceEpoch = newRows[newRows.length - 1]!.epoch + 1;
            const log = await this.api.tool<{ text: string; next_offset: number }>('train_log', {
                job_id: this.jobId,
                offset: this.logOffset,
            });
            let newLog = '';
            if (!isToolError(log)) {
                newLog = log.text;
                this.logOffset = log.next_offset;
            }
            const terminal = isTerminal(status.status);
            let checkpoints: CheckpointInfo[] | null = null;
            if (terminal || newRows.length) {
                const listed = await this.api.tool<{ checkpoints: CheckpointInfo[] }>(
                    'train_checkpoints',
                    { job_id: this.jobId },
                );
                if (!isToolError(listed)) checkpoints = listed.checkpoints;
            }
            this.options.onUpdate({ status, newRows, newLog, checkpoints, terminal });
            if (terminal) {
                this.stop();
                return;
            }
        } catch (e) {
            this.options.onError?.(e instanceof Error ? e.message : String(e));
        }
        this.schedule();
    }
}
