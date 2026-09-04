"""Training jobs: a subprocess supervisor over the run-directory layout.

Pure stdlib, Python 3.8, no torch anywhere -- the manager never imports the
trainer, it **spawns** it (``python -m meshioplusplus.physicsnemo.train
--spec <run_dir>/spec.json``): a long-lived server must not pay the framework
import, and a job must be killable on its own. Each job owns a directory
under ``runs_dir`` (``<YYYYmmdd-HHMMSS>-<4 hex>``) whose files split by
writer: the manager owns ``job.json`` (status, pid, exit code, best
checkpoint); the trainer owns ``progress.json``, ``metrics.jsonl``,
``log.txt`` and ``checkpoints/`` (see ``physicsnemo/_train.py``).

Status is derived, never trusted: a live ``Popen`` handle gives a real exit
code; after a server restart only the pid is left, so liveness is
``os.kill(pid, 0)`` (plus the process start time on Linux, against pid
reuse) and a vanished process with no recorded exit is ``failed`` unless the
trainer's own ``progress.json`` says it completed.

The trainer command can be overridden with the ``MESHIOPLUSPLUS_TRAIN_COMMAND``
environment variable (a JSON list) -- how the tests substitute a tiny fake,
and how a deployment runs the trainer from a different interpreter. It is
deliberately **not** a spec key or a tool parameter: an arbitrary command
reachable over HTTP/MCP would be remote code execution by design.

A ``webhook`` URL (the server's ``--webhook``) is POSTed once per job when it
reaches a terminal state, so a finished run reaches a chat channel or a CI
system without a browser tab staying open. It is posted from the one place a
transition is *observed* -- :meth:`JobManager._refresh` -- and the optional
:class:`Watcher` thread only supplies the polling that makes that observation
happen with nobody asking. Like the trainer command it is a **server-side
setting, never a spec key or a tool parameter**: a client-supplied URL the
server then fetches is server-side request forgery by design.
"""

from __future__ import annotations

import json
import os
import secrets
import shutil
import signal
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from typing import Dict, List, Optional

from ..physicsnemo._train import (
    BEST_CHECKPOINT,
    CARD_SUFFIX,
    CHECKPOINT_DIR,
    JOB_FILE,
    LOG_FILE,
    PROGRESS_FILE,
    SPEC_FILE,
    card_path,
    list_checkpoints,
    read_json,
    read_metrics,
    spec_from_dict,
    spec_to_dict,
    write_json_atomic,
)

TRAIN_COMMAND_ENV = "MESHIOPLUSPLUS_TRAIN_COMMAND"
_ERR = "meshio++: train: "
_TERMINAL = ("finished", "failed", "stopped")
#: How long a webhook POST may take before it is abandoned (best-effort).
WEBHOOK_TIMEOUT = 5.0
#: How often the optional watcher re-derives the state of live jobs.
WATCH_INTERVAL = 5.0


def trainer_command() -> List[str]:
    """The command a job runs, before ``--spec <path>``."""
    override = os.environ.get(TRAIN_COMMAND_ENV)
    if override:
        try:
            command = json.loads(override)
        except ValueError:
            raise ValueError(f"{_ERR}{TRAIN_COMMAND_ENV} must be a JSON list") from None
        if not isinstance(command, list) or not all(
            isinstance(c, str) for c in command
        ):
            raise ValueError(
                f"{_ERR}{TRAIN_COMMAND_ENV} must be a JSON list of strings"
            )
        return command
    return [sys.executable, "-m", "meshioplusplus.physicsnemo.train"]


def new_job_id() -> str:
    return time.strftime("%Y%m%d-%H%M%S") + "-" + secrets.token_hex(2)


def _pid_start(pid: int) -> Optional[str]:
    """Linux: the process start time (clock ticks) from /proc, to tell a
    reused pid from the original; None elsewhere."""
    try:
        with open(f"/proc/{pid}/stat", encoding="utf-8") as fh:
            stat = fh.read()
        # field 22 (1-based) is starttime; the comm field may contain spaces
        return stat[stat.rindex(")") + 2 :].split()[19]
    except (OSError, ValueError, IndexError):
        return None


def is_alive(pid: int, pid_start: Optional[str] = None) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    if pid_start is not None:
        current = _pid_start(pid)
        if current is not None and current != pid_start:
            return False
    return True


def post_webhook(url: str, payload: dict, *, timeout: float = WEBHOOK_TIMEOUT) -> bool:
    """POST ``payload`` as JSON. Best-effort: a failure is logged and dropped,
    never raised into the caller (a job's outcome does not depend on whoever
    wanted to hear about it)."""
    body = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(
        url, data=body, headers={"Content-Type": "application/json"}, method="POST"
    )
    try:
        with urllib.request.urlopen(
            request, timeout=timeout
        ):  # noqa: S310 - operator-supplied
            return True
    except (urllib.error.URLError, OSError, ValueError) as e:
        print(f"{_ERR}webhook POST to {url} failed: {e}", file=sys.stderr)
        return False


class JobManager:
    """Start, watch and stop training jobs under ``runs_dir``."""

    def __init__(self, runs_dir: str, *, webhook: Optional[str] = None) -> None:
        self.runs_dir = os.path.abspath(str(runs_dir))
        self.webhook = webhook
        self._procs: Dict[str, subprocess.Popen] = {}

    # -- paths ------------------------------------------------------------- #
    def job_dir(self, job_id: str) -> str:
        if not job_id or os.sep in job_id or "/" in job_id or job_id in (".", ".."):
            raise ValueError(f"{_ERR}invalid job id '{job_id}'")
        return os.path.join(self.runs_dir, job_id)

    def _job_path(self, job_id: str) -> str:
        return os.path.join(self.job_dir(job_id), JOB_FILE)

    def _read_job(self, job_id: str) -> dict:
        job = read_json(self._job_path(job_id))
        if not isinstance(job, dict):
            raise ValueError(f"{_ERR}no job '{job_id}' under {self.runs_dir}")
        return job

    def _write_job(self, job: dict) -> None:
        write_json_atomic(self._job_path(job["job_id"]), job)

    # -- lifecycle --------------------------------------------------------- #
    def start(self, spec_doc: dict, *, command: Optional[List[str]] = None) -> dict:
        """Validate the spec, lay out the run directory, spawn the trainer."""
        spec = spec_from_dict(spec_doc)  # strict, before anything is created
        job_id = new_job_id()
        run_dir = self.job_dir(job_id)
        while os.path.exists(run_dir):  # pragma: no cover - same second, same hex
            job_id = new_job_id()
            run_dir = self.job_dir(job_id)
        os.makedirs(os.path.join(run_dir, CHECKPOINT_DIR), exist_ok=True)
        doc = spec_to_dict(spec)
        doc["RunDir"] = run_dir
        if not os.path.isabs(doc["Manifest"]):
            doc["Manifest"] = os.path.abspath(doc["Manifest"])
        spec_path = os.path.join(run_dir, SPEC_FILE)
        write_json_atomic(spec_path, doc)
        argv = list(command or trainer_command()) + ["--spec", spec_path]
        log = open(os.path.join(run_dir, LOG_FILE), "ab")
        popen_kwargs = {}
        if os.name == "nt":  # pragma: no cover - windows
            popen_kwargs["creationflags"] = getattr(
                subprocess, "CREATE_NEW_PROCESS_GROUP", 0
            )
        else:
            popen_kwargs["start_new_session"] = True
        try:
            proc = subprocess.Popen(
                argv,
                cwd=run_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL,
                **popen_kwargs,
            )
        finally:
            log.close()
        self._procs[job_id] = proc
        job = {
            "job_id": job_id,
            "status": "running",
            "pid": proc.pid,
            "pid_start": _pid_start(proc.pid),
            "started": time.time(),
            "finished": None,
            "exit_code": None,
            "stop_requested": False,
            "best_checkpoint": None,
            "spec_path": spec_path,
            "manifest": doc["Manifest"],
            "command": argv,
        }
        self._write_job(job)
        return self.status(job_id)

    def _refresh(self, job: dict) -> dict:
        """Derive the current status; persist a transition when one happened."""
        if job["status"] in _TERMINAL:
            return job
        proc = self._procs.get(job["job_id"])
        exit_code = None
        alive = True
        if proc is not None:
            exit_code = proc.poll()
            alive = exit_code is None
        else:
            alive = is_alive(int(job.get("pid") or 0), job.get("pid_start"))
        if alive:
            return job
        progress = read_json(
            os.path.join(self.job_dir(job["job_id"]), PROGRESS_FILE), {}
        )
        completed = bool(isinstance(progress, dict) and progress.get("completed"))
        if job.get("stop_requested"):
            status = "stopped"
        elif exit_code == 0 or (exit_code is None and completed):
            status = "finished"
        else:
            status = "failed"
        job["status"] = status
        job["exit_code"] = exit_code
        job["finished"] = job.get("finished") or time.time()
        if exit_code is None:
            job["reason"] = "process gone" if not completed else "completed"
        self._write_job(job)
        self._procs.pop(job["job_id"], None)
        # The one place a terminal transition is observed, so the one place
        # the webhook fires -- exactly once per job, whoever noticed.
        self._notify(job, progress if isinstance(progress, dict) else {})
        return job

    def _notify(self, job: dict, progress: dict) -> None:
        if not self.webhook:
            return
        payload = {
            "event": f"run.{job['status']}",
            "job_id": job["job_id"],
            "run_dir": self.job_dir(job["job_id"]),
            "status": job["status"],
            "exit_code": job.get("exit_code"),
            "reason": job.get("reason"),
            "manifest": job.get("manifest"),
            "epoch": progress.get("epoch"),
            "epochs": progress.get("epochs"),
            "best_valid_loss": progress.get("best_valid_loss"),
            "best_checkpoint": job.get("best_checkpoint")
            or progress.get("best_checkpoint"),
        }
        url = self.webhook
        # Off the caller's thread: a slow endpoint must not stall the request
        # (or the watcher) that happened to notice the transition.
        threading.Thread(
            target=post_webhook,
            args=(url, payload),
            daemon=True,
            name="meshioplusplus-webhook",
        ).start()

    def status(self, job_id: str) -> dict:
        job = self._refresh(self._read_job(job_id))
        run_dir = self.job_dir(job_id)
        progress = read_json(os.path.join(run_dir, PROGRESS_FILE), {})
        if not isinstance(progress, dict):
            progress = {}
        rows = read_metrics(run_dir)
        return {
            "job_id": job_id,
            "run_dir": run_dir,
            "status": job["status"],
            "pid": job.get("pid"),
            "started": job.get("started"),
            "finished": job.get("finished"),
            "exit_code": job.get("exit_code"),
            "reason": job.get("reason"),
            "manifest": job.get("manifest"),
            "best_checkpoint": job.get("best_checkpoint")
            or progress.get("best_checkpoint"),
            "epoch": progress.get("epoch"),
            "epochs": progress.get("epochs"),
            "best_epoch": progress.get("best_epoch"),
            "best_valid_loss": progress.get("best_valid_loss"),
            "eta_seconds": progress.get("eta_seconds"),
            "device": progress.get("device"),
            "completed": bool(progress.get("completed")),
            "num_metrics": len(rows),
            "last": rows[-1] if rows else None,
        }

    def list_jobs(
        self, status: Optional[str] = None, manifest: Optional[str] = None
    ) -> List[dict]:
        """Every job under ``runs_dir`` (newest first), each summarized."""
        out = []
        if not os.path.isdir(self.runs_dir):
            return out
        for name in sorted(os.listdir(self.runs_dir), reverse=True):
            if not os.path.isfile(os.path.join(self.runs_dir, name, JOB_FILE)):
                continue
            try:
                summary = self._summary(name)
            except ValueError:
                continue
            if status and summary["status"] != status:
                continue
            if manifest and os.path.abspath(
                summary.get("manifest") or ""
            ) != os.path.abspath(manifest):
                continue
            out.append(summary)
        return out

    def _summary(self, job_id: str) -> dict:
        state = self.status(job_id)
        spec = read_json(os.path.join(self.job_dir(job_id), SPEC_FILE), {}) or {}
        rows = read_metrics(self.job_dir(job_id))
        last = rows[-1] if rows else None
        finite = [r["valid_loss"] for r in rows if r.get("valid_loss") is not None]
        state.update(
            {
                "fields": spec.get("Fields", []),
                "target_fields": spec.get("TargetFields", []),
                "train_split": spec.get("TrainSplit"),
                "valid_split": spec.get("ValidSplit"),
                "epochs": spec.get("Epochs", state.get("epochs")),
                "batch_size": spec.get("BatchSize"),
                "learning_rate": spec.get("LearningRate"),
                "seed": spec.get("Seed"),
                "hidden_dim": (spec.get("Model") or {}).get("HiddenDim"),
                "processor_size": (spec.get("Model") or {}).get("ProcessorSize"),
                "tags": spec.get("Tags", []),
                "notes": spec.get("Notes"),
                "final_train_loss": last["train_loss"] if last else None,
                "final_valid_loss": last["valid_loss"] if last else None,
                "best_valid_loss": (
                    min(finite) if finite else state.get("best_valid_loss")
                ),
                "duration_seconds": (
                    (state["finished"] or time.time()) - state["started"]
                    if state.get("started")
                    else None
                ),
            }
        )
        return state

    def stop(self, job_id: str, grace: float = 10.0) -> dict:
        """SIGTERM the job's process group (the trainer finishes its epoch
        and writes ``final.*``), SIGKILL after ``grace`` seconds."""
        job = self._refresh(self._read_job(job_id))
        if job["status"] in _TERMINAL:
            return self.status(job_id)
        job["stop_requested"] = True
        self._write_job(job)
        pid = int(job.get("pid") or 0)
        proc = self._procs.get(job_id)
        try:
            if os.name == "nt":  # pragma: no cover - windows
                if proc is not None:
                    proc.send_signal(
                        getattr(signal, "CTRL_BREAK_EVENT", signal.SIGTERM)
                    )
                else:
                    os.kill(pid, signal.SIGTERM)
            else:
                os.killpg(pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError, OSError):
            pass
        deadline = time.time() + max(float(grace), 0.0)
        while time.time() < deadline:
            if proc is not None:
                if proc.poll() is not None:
                    break
            elif not is_alive(pid, job.get("pid_start")):
                break
            time.sleep(0.05)
        else:
            try:
                if os.name == "nt":  # pragma: no cover - windows
                    if proc is not None:
                        proc.kill()
                    else:
                        os.kill(pid, signal.SIGTERM)
                else:
                    os.killpg(pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError, OSError):
                pass
            if proc is not None:
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:  # pragma: no cover
                    pass
        return self.status(job_id)

    # -- files ------------------------------------------------------------- #
    def log(self, job_id: str, offset: int = 0, max_bytes: int = 65536) -> dict:
        """A window of ``log.txt`` from ``offset``; ``done`` once the job is
        terminal and the window reached the end."""
        state = self.status(job_id)
        path = os.path.join(self.job_dir(job_id), LOG_FILE)
        offset = max(int(offset), 0)
        chunk = b""
        size = 0
        try:
            size = os.path.getsize(path)
            with open(path, "rb") as fh:
                fh.seek(offset)
                chunk = fh.read(max(int(max_bytes), 0))
        except OSError:
            pass
        next_offset = offset + len(chunk)
        return {
            "job_id": job_id,
            "status": state["status"],
            "text": chunk.decode("utf-8", errors="replace"),
            "offset": offset,
            "next_offset": next_offset,
            "size": size,
            "done": state["status"] in _TERMINAL and next_offset >= size,
        }

    def metrics(self, job_id: str, since_epoch: int = 0) -> dict:
        state = self.status(job_id)
        return {
            "job_id": job_id,
            "status": state["status"],
            "rows": read_metrics(self.job_dir(job_id), int(since_epoch)),
        }

    def checkpoints(self, job_id: str) -> dict:
        state = self.status(job_id)
        return {
            "job_id": job_id,
            "status": state["status"],
            "best_checkpoint": state["best_checkpoint"],
            "checkpoints": list_checkpoints(
                self.job_dir(job_id), state["best_checkpoint"]
            ),
        }

    def mark_best(self, job_id: str, checkpoint: str) -> dict:
        """Copy a checkpoint (and its card) to ``best.mdlus`` and record it."""
        job = self._read_job(job_id)
        directory = os.path.join(self.job_dir(job_id), CHECKPOINT_DIR)
        source = os.path.join(directory, os.path.basename(str(checkpoint)))
        if not os.path.isfile(source) or not source.endswith(".mdlus"):
            raise ValueError(f"{_ERR}no checkpoint '{checkpoint}' in job '{job_id}'")
        target = os.path.join(directory, BEST_CHECKPOINT)
        if os.path.abspath(source) != os.path.abspath(target):
            shutil.copyfile(source, target)
            if os.path.isfile(card_path(source)):
                shutil.copyfile(card_path(source), card_path(target))
            else:
                try:
                    os.remove(card_path(target))
                except OSError:
                    pass
        job["best_checkpoint"] = target
        job["best_source"] = os.path.basename(source)
        self._write_job(job)
        return self.checkpoints(job_id)

    def resolve_checkpoint(self, job_id: str, checkpoint: Optional[str] = None) -> str:
        """The checkpoint a job should predict with: an explicit name, else
        the marked best, else ``best.mdlus``, else ``final.mdlus``."""
        directory = os.path.join(self.job_dir(job_id), CHECKPOINT_DIR)
        if checkpoint:
            path = os.path.join(directory, os.path.basename(str(checkpoint)))
            if os.path.isfile(path):
                return path
            raise ValueError(f"{_ERR}no checkpoint '{checkpoint}' in job '{job_id}'")
        job = self._read_job(job_id)
        for candidate in (
            job.get("best_checkpoint"),
            os.path.join(directory, BEST_CHECKPOINT),
            os.path.join(directory, "final.mdlus"),
        ):
            if candidate and os.path.isfile(candidate):
                return candidate
        raise ValueError(f"{_ERR}job '{job_id}' has no checkpoint yet")


class Watcher:
    """Poll live jobs so a terminal transition is noticed (and its webhook
    posted) with nobody asking -- a run killed outright still notifies.

    A daemon thread calling :meth:`JobManager.status`, which is the same
    derivation every request performs; it adds no state of its own.
    """

    def __init__(
        self, manager: JobManager, *, interval: float = WATCH_INTERVAL
    ) -> None:
        self.manager = manager
        self.interval = float(interval)
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None

    def start(self) -> "Watcher":
        if self._thread is not None:
            return self
        self._thread = threading.Thread(
            target=self._run, daemon=True, name="meshioplusplus-job-watcher"
        )
        self._thread.start()
        return self

    def stop(self, timeout: Optional[float] = None) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout)
            self._thread = None

    def poll_once(self) -> int:
        """Re-derive every non-terminal job; returns how many were checked."""
        checked = 0
        try:
            names = sorted(os.listdir(self.manager.runs_dir))
        except OSError:
            return 0
        for name in names:
            if not os.path.isfile(os.path.join(self.manager.runs_dir, name, JOB_FILE)):
                continue
            try:
                job = self.manager._read_job(name)
                if job.get("status") in _TERMINAL:
                    continue
                self.manager.status(name)
                checked += 1
            except Exception as e:  # noqa: BLE001 - a bad job must not stop the watch
                print(f"{_ERR}watcher: {name}: {e}", file=sys.stderr)
        return checked

    def _run(self) -> None:
        while not self._stop.wait(self.interval):
            self.poll_once()


__all__ = [
    "CARD_SUFFIX",
    "JobManager",
    "TRAIN_COMMAND_ENV",
    "WATCH_INTERVAL",
    "WEBHOOK_TIMEOUT",
    "Watcher",
    "is_alive",
    "new_job_id",
    "post_webhook",
    "trainer_command",
]
