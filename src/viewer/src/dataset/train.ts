/**
 * Launching and following a training run from the dataset manager.
 *
 * Everything here needs the [companion process](doc/dashboard.md): real GPU
 * training cannot happen in a browser tab, so the page's job is to describe
 * the run, start it (`train_start`), and then *follow* it — the poller in
 * `jobs.ts` over the incremental `train_status`/`train_metrics`/`train_log`
 * tools. The panel owns no training logic of its own; the spec it builds is
 * the same PascalCase document `python -m meshioplusplus.physicsnemo.train`
 * takes, and the server writes it into the run directory.
 */

import { isToolError, type ApiClient } from './api';
import { SERIES_COLORS, renderLineChart, type ChartSeries } from './chart';
import { JobPoller, isTerminal, type CheckpointInfo, type JobStatus, type JobSummary, type MetricRow } from './jobs';
import { notifyRunFinished, requestNotifications, resetTitle } from './notify';
import { $, setOptions, show, svg$ } from '../ui/dom';
import type { ActiveJob, DatasetState } from './types';

/** Bytes -> a short human size for the checkpoint list. */
function humanSize(bytes: number): string {
    const units = ['B', 'kB', 'MB', 'GB'];
    let value = bytes;
    let unit = 0;
    while (value >= 1024 && unit < units.length - 1) {
        value /= 1024;
        unit += 1;
    }
    return `${value < 10 && unit ? value.toFixed(1) : Math.round(value)} ${units[unit]}`;
}

function humanSeconds(seconds: number | null): string {
    if (seconds === null || !Number.isFinite(seconds)) return '—';
    const total = Math.max(Math.round(seconds), 0);
    const h = Math.floor(total / 3600);
    const m = Math.floor((total % 3600) / 60);
    const s = total % 60;
    if (h) return `${h}h ${m}m`;
    if (m) return `${m}m ${s}s`;
    return `${s}s`;
}

function multiValues(select: HTMLSelectElement): string[] {
    return [...select.selectedOptions].map((o) => o.value);
}

function setMulti(select: HTMLSelectElement, names: string[], selected: string[]): void {
    select.replaceChildren(
        ...names.map((name) => {
            const option = document.createElement('option');
            option.value = name;
            option.textContent = name;
            option.selected = selected.includes(name);
            return option;
        }),
    );
    select.size = Math.min(Math.max(names.length, 3), 6);
}

export interface TrainPanelHooks {
    /** The connected client, or null. */
    getApi(): ApiClient | null;
    /** The open manifest's path ON THE SERVER (bound by content hash), or
     * null when no manifest is open or it is not bound. */
    getManifestPath(): string | null;
    /** The open manifest's entries, for the prediction picker. */
    getEntries(): { id: string; split: string | null }[];
    /** Render a mesh the server produced (a prediction), colouring it by its
     * error array — the viewer's own session slot, so the staged dataset
     * entry is untouched. */
    previewPrediction(name: string, bytes: ArrayBuffer): Promise<void>;
    setState(patch: Partial<DatasetState>): void;
    setStatus(message: string): void;
    fail(error: unknown): void;
}

export class TrainPanel {
    private poller: JobPoller | null = null;
    private active: ActiveJob | null = null;
    private jobs: JobSummary[] = [];
    private fieldNames: string[] = [];
    private logText = '';

    constructor(private readonly hooks: TrainPanelHooks) {
        $('t-start').addEventListener('click', () => void this.start().catch(hooks.fail));
        $('t-model').addEventListener('change', () => this.syncModelOptions());
        this.syncModelOptions();
        $('run-close').addEventListener('click', () => this.close());
        $('pred-run').addEventListener('click', () => void this.predict().catch(hooks.fail));
        $('run-stop').addEventListener('click', () => void this.stop().catch(hooks.fail));
        $<HTMLSelectElement>('run-job').addEventListener('change', (e) => {
            const jobId = (e.target as HTMLSelectElement).value;
            if (jobId) void this.follow(jobId).catch(hooks.fail);
        });
    }

    /** Enable/disable the form for the manifest now open. */
    refreshAvailability(): void {
        const api = this.hooks.getApi();
        const manifest = this.hooks.getManifestPath();
        const ready = !!api && !!manifest;
        $<HTMLButtonElement>('t-start').disabled = !ready;
        $<HTMLButtonElement>('t-runs').disabled = !api;
        $('t-hint').textContent = api
            ? ready
                ? 'Trains on the companion process; follow it here.'
                : 'Save the manifest so the server can see it, then reload the datasets.'
            : 'Connect a companion process to train.';
        if (ready) void this.refreshDefaults().catch(this.hooks.fail);
    }

    /** Populate the field/target/split pickers from the manifest. */
    async refreshDefaults(): Promise<void> {
        const api = this.hooks.getApi();
        const manifest = this.hooks.getManifestPath();
        if (!api || !manifest) return;
        const report = await api.tool<{
            available_fields: { point: string[]; cell: string[] };
            splits: Record<string, number>;
            spec: Record<string, unknown>;
            frameworks: Record<string, boolean>;
        }>('train_defaults', { manifest_path: manifest });
        if (isToolError(report)) throw new Error(report.error);
        this.fieldNames = report.available_fields.point;
        const previousFields = multiValues($<HTMLSelectElement>('t-fields'));
        const previousTargets = multiValues($<HTMLSelectElement>('t-targets'));
        setMulti($<HTMLSelectElement>('t-fields'), this.fieldNames, previousFields);
        setMulti($<HTMLSelectElement>('t-targets'), this.fieldNames, previousTargets);
        const splits = Object.keys(report.splits).filter(Boolean).sort();
        setOptions(
            $<HTMLSelectElement>('t-train-split'),
            splits.map((s) => ({ value: s, label: s })),
            splits.includes('train') ? 'train' : (splits[0] ?? ''),
        );
        setOptions(
            $<HTMLSelectElement>('t-valid-split'),
            splits.map((s) => ({ value: s, label: s })),
            splits.includes('valid') ? 'valid' : (splits[1] ?? splits[0] ?? ''),
        );
        const missing = Object.entries(report.frameworks)
            .filter(([, present]) => !present)
            .map(([name]) => name);
        if (missing.length) {
            $('t-error').textContent =
                `the companion process has no ${missing.join(' / ')} installed; ` +
                'starting a run will fail by name until it does';
            $('t-error').hidden = false;
        } else {
            $('t-error').hidden = true;
        }
    }

    /** The jobs the server knows; `open` also shows the newest one. */
    async refreshJobs(open = false): Promise<void> {
        const api = this.hooks.getApi();
        if (!api) return;
        const listed = await api.tool<{ jobs: JobSummary[] }>('train_list', {});
        if (isToolError(listed)) throw new Error(listed.error);
        this.jobs = listed.jobs;
        this.hooks.setState({ jobs: this.jobs });
        setOptions(
            $<HTMLSelectElement>('run-job'),
            this.jobs.map((job) => ({
                value: job.job_id,
                label: `${job.job_id} · ${job.status}`,
            })),
            this.active?.jobId ?? this.jobs[0]?.job_id,
        );
        if (open) {
            const jobId = this.active?.jobId ?? this.jobs[0]?.job_id;
            if (jobId) await this.follow(jobId);
            else this.hooks.setStatus('no training runs yet');
        }
    }

    private specFromForm(manifest: string): Record<string, unknown> {
        const fields = multiValues($<HTMLSelectElement>('t-fields'));
        const targets = multiValues($<HTMLSelectElement>('t-targets'));
        if (!fields.length || !targets.length) {
            throw new Error('meshio++: pick at least one input field and one target field');
        }
        const number = (id: string) => Number($<HTMLInputElement>(id).value);
        const modelName = $<HTMLSelectElement>('t-model').value;
        const common = {
            manifest_path: manifest,
            fields,
            target_fields: targets,
            train_split: $<HTMLSelectElement>('t-train-split').value,
            valid_split: $<HTMLSelectElement>('t-valid-split').value,
            epochs: number('t-epochs'),
            batch_size: number('t-batch'),
            learning_rate: number('t-lr'),
            seed: number('t-seed'),
            model_name: modelName,
        };
        // Only the chosen family's hyperparameters are sent. The server refuses
        // the other family's by name, so posting both would turn a UI default
        // into an error the user never asked for.
        if (modelName === 'srresnet') {
            const resolution = $<HTMLInputElement>('t-resolution')
                .value.split(',')
                .map((part) => Number(part.trim()));
            if (resolution.length !== 3 || resolution.some((v) => !(v > 0))) {
                throw new Error('meshio++: resolution must be three positive cell counts, e.g. 16,16,16');
            }
            return {
                ...common,
                resolution,
                scaling_factor: Number($<HTMLSelectElement>('t-scaling').value),
                conv_layer_size: number('t-conv'),
                resid_blocks: number('t-blocks'),
            };
        }
        return {
            ...common,
            processor_size: number('t-processor'),
            hidden_dim: number('t-hidden'),
        };
    }

    /** Show only the chosen family's options. */
    syncModelOptions(): void {
        const isGrid = $<HTMLSelectElement>('t-model').value === 'srresnet';
        $('t-graph-opts').hidden = isGrid;
        $('t-grid-opts').hidden = !isGrid;
    }

    async start(): Promise<void> {
        const api = this.hooks.getApi();
        const manifest = this.hooks.getManifestPath();
        if (!api || !manifest) return;
        $('t-error').hidden = true;
        let args: Record<string, unknown>;
        try {
            args = this.specFromForm(manifest);
        } catch (e) {
            $('t-error').textContent = e instanceof Error ? e.message : String(e);
            $('t-error').hidden = false;
            return;
        }
        // The permission prompt needs this user gesture, but the run must not
        // wait on the answer: a prompt left open would otherwise block
        // starting a run, and the answer only decides whether the finish is
        // announced.
        void requestNotifications().then((permission) =>
            this.hooks.setState({ notifications: permission }),
        );
        this.hooks.setStatus('starting the run…');
        const started = await api.tool<JobStatus>('train_start', args);
        if (isToolError(started)) {
            $('t-error').textContent = started.error;
            $('t-error').hidden = false;
            this.hooks.setStatus('ready');
            return;
        }
        this.hooks.setStatus(`run ${started.job_id} started`);
        await this.refreshJobs();
        await this.follow(started.job_id);
    }

    /** Predict one entry with this run's checkpoint and show the result in
     * the mesh viewer, coloured by its own error field. */
    async predict(): Promise<void> {
        const api = this.hooks.getApi();
        const manifest = this.hooks.getManifestPath();
        const job = this.active;
        const entryId = $<HTMLSelectElement>('pred-entry').value;
        if (!api || !manifest || !job || !entryId) return;
        const metrics = $('pred-metrics');
        metrics.textContent = 'predicting…';
        this.hooks.setStatus(`predicting ${entryId}…`);
        const report = await api.tool<{
            checkpoint: string;
            predictions: { entry_id: string; output_path: string; rmse: number | null; max_error: number | null }[];
        }>('train_predict', { manifest_path: manifest, job_id: job.jobId, split: null, entry_ids: [entryId] });
        if (isToolError(report)) {
            metrics.textContent = report.error;
            this.hooks.setStatus('ready');
            return;
        }
        const row = report.predictions[0];
        if (!row) {
            metrics.textContent = 'the server predicted nothing for that entry';
            return;
        }
        // The prediction is a file on the server: fetch it through the same
        // sandboxed route the checkpoint downloads use, then render it.
        const response = await fetch(api.fileUrl(row.output_path));
        if (!response.ok) throw new Error(`meshio++: could not fetch ${row.output_path}`);
        await this.hooks.previewPrediction(`${row.entry_id}.vtu`, await response.arrayBuffer());
        const rmse = row.rmse === null ? '—' : row.rmse.toPrecision(3);
        const max = row.max_error === null ? '—' : row.max_error.toPrecision(3);
        metrics.textContent = `${row.entry_id}: RMSE ${rmse} · max |error| ${max}`;
        this.hooks.setState({
            prediction: {
                jobId: job.jobId,
                entryId: row.entry_id,
                outputPath: row.output_path,
                rmse: row.rmse,
                maxError: row.max_error,
            },
        });
        this.hooks.setStatus(`predicted ${row.entry_id}`);
    }

    /** Fill the prediction picker from the open manifest. */
    private refreshEntries(): void {
        const entries = this.hooks.getEntries();
        setOptions(
            $<HTMLSelectElement>('pred-entry'),
            entries.map((e) => ({
                value: e.id,
                label: e.split ? `${e.id} (${e.split})` : e.id,
            })),
            // Default to a held-out entry: predicting on a training case is
            // rarely what you want to look at.
            entries.find((e) => e.split && e.split !== 'train')?.id ?? entries[0]?.id,
        );
        $<HTMLButtonElement>('pred-run').disabled = entries.length === 0;
    }

    /** Show a job's panel and poll it while it runs. */
    async follow(jobId: string): Promise<void> {
        const api = this.hooks.getApi();
        if (!api) return;
        this.poller?.stop();
        this.logText = '';
        $('run-log').textContent = '';
        this.active = {
            jobId,
            status: 'running',
            epoch: 0,
            epochs: 0,
            metrics: [],
            checkpoints: [],
            bestCheckpoint: null,
            error: null,
        };
        $<HTMLSelectElement>('run-job').value = jobId;
        $('pred-metrics').textContent = '';
        this.refreshEntries();
        show($('run-wrap'));
        this.render();
        this.poller = new JobPoller(api, jobId, {
            onUpdate: (update) => this.onUpdate(update.status, update.newRows, update.newLog, update.checkpoints, update.terminal),
            onError: (message) => {
                if (this.active) this.active.error = message;
                this.render();
            },
        });
        this.poller.start();
    }

    private onUpdate(
        status: JobStatus,
        newRows: MetricRow[],
        newLog: string,
        checkpoints: CheckpointInfo[] | null,
        terminal: boolean,
    ): void {
        if (!this.active || this.active.jobId !== status.job_id) return;
        this.active.status = status.status;
        this.active.epoch = status.epoch ?? this.active.epoch;
        this.active.epochs = status.epochs ?? this.active.epochs;
        this.active.bestCheckpoint = status.best_checkpoint;
        this.active.error = null;
        this.active.metrics = [...this.active.metrics, ...newRows];
        if (checkpoints) this.active.checkpoints = checkpoints;
        if (newLog) {
            this.logText += newLog;
            const pre = $('run-log');
            pre.textContent = this.logText;
            if ($<HTMLInputElement>('run-follow').checked) pre.scrollTop = pre.scrollHeight;
        }
        this.render(status);
        if (terminal) {
            const detail = status.best_valid_loss !== null
                ? `best validation loss ${status.best_valid_loss.toPrecision(3)}`
                : (status.reason ?? '');
            notifyRunFinished(status.job_id, status.status, detail);
            this.hooks.setStatus(`run ${status.job_id} ${status.status}`);
            void this.refreshJobs().catch(() => undefined);
        }
    }

    async stop(): Promise<void> {
        const api = this.hooks.getApi();
        if (!api || !this.active) return;
        this.hooks.setStatus(`stopping ${this.active.jobId}…`);
        const stopped = await api.tool<JobStatus>('train_stop', { job_id: this.active.jobId });
        if (isToolError(stopped)) throw new Error(stopped.error);
        this.onUpdate(stopped, [], '', null, isTerminal(stopped.status));
    }

    close(): void {
        this.poller?.stop();
        this.poller = null;
        $('run-wrap').hidden = true;
        resetTitle();
    }

    dispose(): void {
        this.close();
        this.active = null;
        this.jobs = [];
        this.hooks.setState({ jobs: [], activeJob: null, compare: [], prediction: null });
    }

    private render(status?: JobStatus): void {
        const job = this.active;
        if (!job) return;
        $('run-title').textContent = job.jobId;
        const state = status?.status ?? job.status;
        const badge = $('run-status');
        badge.textContent = state + (status?.reason ? ` (${status.reason})` : '');
        badge.className = `badge ${state === 'finished' ? 'ok' : state === 'running' ? 'none' : 'bad'}`;
        const progress = $<HTMLProgressElement>('run-progress');
        progress.max = Math.max(job.epochs, 1);
        progress.value = job.epoch;
        const parts = [`epoch ${job.epoch}/${job.epochs || '?'}`];
        if (status?.eta_seconds !== null && status?.eta_seconds !== undefined && !isTerminal(state)) {
            parts.push(`ETA ${humanSeconds(status.eta_seconds)}`);
        }
        if (status?.best_valid_loss !== null && status?.best_valid_loss !== undefined) {
            parts.push(`best ${status.best_valid_loss.toPrecision(3)} @ ${status.best_epoch}`);
        }
        if (status?.device) parts.push(status.device);
        if (job.error) parts.push(job.error);
        $('run-eta').textContent = parts.join(' · ');
        $<HTMLButtonElement>('run-stop').disabled = isTerminal(state);

        const series: ChartSeries[] = [
            {
                name: 'train',
                color: SERIES_COLORS[0]!,
                points: job.metrics.map((r) => [r.epoch, r.train_loss] as [number, number]),
            },
        ];
        const validPoints = job.metrics
            .filter((r) => r.valid_loss !== null)
            .map((r) => [r.epoch, r.valid_loss as number] as [number, number]);
        if (validPoints.length) {
            series.push({ name: 'valid', color: SERIES_COLORS[1]!, points: validPoints });
        }
        renderLineChart(svg$('run-chart'), series, {
            logY: true,
            xLabel: 'epoch',
            yLabel: 'loss (normalized MSE)',
            tooltip: $('run-tooltip'),
            formatX: (x) => `epoch ${Math.round(x)}`,
        });
        $('run-legend').replaceChildren(
            ...series.map((s) => {
                const item = document.createElement('span');
                item.className = 'split-item';
                const swatch = document.createElement('i');
                swatch.className = 'split-swatch';
                swatch.style.background = s.color;
                item.append(swatch, s.name);
                return item;
            }),
        );
        this.renderCheckpoints();
        this.hooks.setState({ activeJob: { ...job, metrics: [...job.metrics] } });
    }

    private renderCheckpoints(): void {
        const api = this.hooks.getApi();
        const job = this.active;
        const list = $('run-checkpoints');
        if (!job || !api) {
            list.replaceChildren();
            return;
        }
        list.replaceChildren(
            ...job.checkpoints.map((checkpoint) => {
                const li = document.createElement('li');
                if (checkpoint.is_best) li.classList.add('best');
                const label = document.createElement('span');
                label.className = 'ckpt-name';
                const loss =
                    checkpoint.valid_loss === null ? '' : ` · ${checkpoint.valid_loss.toPrecision(3)}`;
                label.textContent =
                    `${checkpoint.name} · ${checkpoint.kind}` +
                    (checkpoint.epoch === null ? '' : ` · epoch ${checkpoint.epoch}`) +
                    `${loss} · ${humanSize(checkpoint.size)}`;
                const download = document.createElement('a');
                download.href = api.fileUrl(checkpoint.path);
                download.textContent = 'download';
                download.setAttribute('download', checkpoint.name);
                li.append(label, download);
                if (!checkpoint.is_best) {
                    const mark = document.createElement('button');
                    mark.type = 'button';
                    mark.textContent = 'mark best';
                    mark.addEventListener('click', () => void this.markBest(checkpoint.name));
                    li.append(mark);
                }
                return li;
            }),
        );
    }

    private async markBest(name: string): Promise<void> {
        const api = this.hooks.getApi();
        const job = this.active;
        if (!api || !job) return;
        const marked = await api.tool<{ checkpoints: CheckpointInfo[]; best_checkpoint: string }>(
            'train_mark_best',
            { job_id: job.jobId, checkpoint: name },
        );
        if (isToolError(marked)) {
            this.hooks.fail(new Error(marked.error));
            return;
        }
        job.checkpoints = marked.checkpoints;
        job.bestCheckpoint = marked.best_checkpoint;
        this.render();
        this.hooks.setStatus(`marked ${name} best`);
    }
}
