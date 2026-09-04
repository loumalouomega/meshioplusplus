"""The training job manager (``meshioplusplus.mcp._jobs``) over the fake
trainer in ``tests/python/_fake_trainer.py`` -- pure stdlib, no torch, runs
in the default matrix. Signals are POSIX; the stop tests skip on Windows.
"""

from __future__ import annotations

import json
import os
import sys
import time

import pytest

from meshioplusplus.mcp import _jobs
from meshioplusplus.physicsnemo._train import (
    JOB_FILE,
    PROGRESS_FILE,
    read_json,
    write_json_atomic,
)

FAKE = [sys.executable, os.path.join(os.path.dirname(__file__), "_fake_trainer.py")]
posix_only = pytest.mark.skipif(os.name == "nt", reason="process groups are POSIX here")


@pytest.fixture()
def manager(tmp_path, monkeypatch):
    monkeypatch.setenv(_jobs.TRAIN_COMMAND_ENV, json.dumps(FAKE))
    return _jobs.JobManager(str(tmp_path / "runs"))


def _spec(tmp_path, **extra):
    manifest = tmp_path / "m.json"
    manifest.write_text('{"Version": 1, "Entries": []}', encoding="utf-8")
    return {
        "Manifest": str(manifest),
        "Fields": ["q"],
        "TargetFields": ["T"],
        "Epochs": 4,
        "CheckpointEvery": 2,
        **extra,
    }


def _wait(manager, job_id, until=("finished", "failed", "stopped"), timeout=20.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        state = manager.status(job_id)
        if state["status"] in until:
            return state
        time.sleep(0.05)
    raise AssertionError(
        f"job {job_id} never reached {until}: {manager.status(job_id)}"
    )


def test_trainer_command_env_override(monkeypatch):
    monkeypatch.delenv(_jobs.TRAIN_COMMAND_ENV, raising=False)
    assert _jobs.trainer_command() == [
        sys.executable,
        "-m",
        "meshioplusplus.physicsnemo.train",
    ]
    monkeypatch.setenv(_jobs.TRAIN_COMMAND_ENV, '["python", "x.py"]')
    assert _jobs.trainer_command() == ["python", "x.py"]
    monkeypatch.setenv(_jobs.TRAIN_COMMAND_ENV, "not json")
    with pytest.raises(ValueError, match="JSON list"):
        _jobs.trainer_command()
    monkeypatch.setenv(_jobs.TRAIN_COMMAND_ENV, "[1]")
    with pytest.raises(ValueError, match="JSON list of strings"):
        _jobs.trainer_command()


def test_start_runs_to_completion_with_metrics_log_and_checkpoints(manager, tmp_path):
    started = manager.start(_spec(tmp_path, Tags=["t"]))
    job_id = started["job_id"]
    assert started["status"] == "running" and started["pid"] > 0
    assert os.path.isfile(os.path.join(started["run_dir"], "spec.json"))
    state = _wait(manager, job_id)
    assert state["status"] == "finished" and state["exit_code"] == 0
    assert state["completed"] and state["epoch"] == 4 and state["epochs"] == 4
    assert state["num_metrics"] == 4 and state["last"]["epoch"] == 3
    assert state["best_checkpoint"].endswith("best.mdlus")
    log = manager.log(job_id)
    assert "fake trainer: done" in log["text"] and log["done"]
    head = manager.log(job_id, offset=0, max_bytes=10)
    assert len(head["text"]) == 10 and head["next_offset"] == 10 and not head["done"]
    rest = manager.log(job_id, offset=head["next_offset"])
    assert head["text"] + rest["text"] == log["text"]
    rows = manager.metrics(job_id, since_epoch=2)["rows"]
    assert [r["epoch"] for r in rows] == [2, 3]
    ckpts = manager.checkpoints(job_id)
    names = [c["name"] for c in ckpts["checkpoints"]]
    assert (
        "Fake.0.1.mdlus" in names
        and "Fake.0.3.mdlus" in names
        and "final.mdlus" in names
    )
    assert ckpts["best_checkpoint"].endswith("best.mdlus")
    assert [c for c in ckpts["checkpoints"] if c["is_best"]][0]["name"] == "best.mdlus"
    summary = manager.list_jobs()[0]
    assert summary["job_id"] == job_id and summary["tags"] == ["t"]
    assert summary["final_valid_loss"] == pytest.approx(1.2 / 4)
    assert summary["best_valid_loss"] == pytest.approx(1.2 / 4)
    assert manager.list_jobs(status="running") == []
    assert manager.list_jobs(manifest=str(tmp_path / "m.json"))[0]["job_id"] == job_id
    assert manager.list_jobs(manifest="/nope.json") == []


def test_a_failing_trainer_is_failed_with_its_exit_code(manager, tmp_path):
    job_id = manager.start(_spec(tmp_path, Notes="fail"))["job_id"]
    state = _wait(manager, job_id)
    assert state["status"] == "failed" and state["exit_code"] == 3
    assert "failing on purpose" in manager.log(job_id)["text"]


@posix_only
def test_stop_terminates_the_group_and_records_stopped(manager, tmp_path):
    job_id = manager.start(_spec(tmp_path, Epochs=200))["job_id"]
    _wait(manager, job_id, until=("running",))
    while manager.status(job_id)["num_metrics"] < 1:
        time.sleep(0.05)
    state = manager.stop(job_id, grace=10.0)
    assert state["status"] == "stopped"
    assert manager.status(job_id)["exit_code"] == 143
    assert "fake trainer: stopped" in manager.log(job_id)["text"]
    assert "final.mdlus" in [
        c["name"] for c in manager.checkpoints(job_id)["checkpoints"]
    ]
    # stopping again is a no-op
    assert manager.stop(job_id)["status"] == "stopped"


@posix_only
def test_stop_escalates_to_kill_when_sigterm_is_ignored(manager, tmp_path):
    job_id = manager.start(_spec(tmp_path, Epochs=200, Notes="hang"))["job_id"]
    while manager.status(job_id)["num_metrics"] < 1:
        time.sleep(0.05)
    state = manager.stop(job_id, grace=0.5)
    assert state["status"] == "stopped"
    assert manager.status(job_id)["exit_code"] not in (0, None)


def test_orphan_after_a_restart_is_derived_from_the_files(manager, tmp_path):
    job_id = manager.start(_spec(tmp_path))["job_id"]
    _wait(manager, job_id)
    fresh = _jobs.JobManager(manager.runs_dir)  # no Popen handles
    # rewrite as "running" with a dead pid: completed progress -> finished
    job = read_json(os.path.join(fresh.job_dir(job_id), JOB_FILE))
    job.update(
        status="running", exit_code=None, finished=None, pid=2**22 - 1, pid_start=None
    )
    write_json_atomic(os.path.join(fresh.job_dir(job_id), JOB_FILE), job)
    state = fresh.status(job_id)
    assert state["status"] == "finished" and state["reason"] == "completed"
    # ...and with no completion recorded -> failed, "process gone"
    job.update(status="running", exit_code=None, finished=None)
    write_json_atomic(os.path.join(fresh.job_dir(job_id), JOB_FILE), job)
    progress = read_json(os.path.join(fresh.job_dir(job_id), PROGRESS_FILE))
    progress["completed"] = False
    write_json_atomic(os.path.join(fresh.job_dir(job_id), PROGRESS_FILE), progress)
    state = fresh.status(job_id)
    assert state["status"] == "failed" and state["reason"] == "process gone"


def test_mark_best_and_resolve_checkpoint(manager, tmp_path):
    job_id = manager.start(_spec(tmp_path))["job_id"]
    _wait(manager, job_id)
    ckpts = manager.mark_best(job_id, "Fake.0.1.mdlus")
    best = [c for c in ckpts["checkpoints"] if c["is_best"]]
    assert len(best) == 1 and best[0]["name"] == "best.mdlus" and best[0]["epoch"] == 1
    assert manager.resolve_checkpoint(job_id).endswith("best.mdlus")
    assert manager.resolve_checkpoint(job_id, "final.mdlus").endswith("final.mdlus")
    with pytest.raises(ValueError, match="no checkpoint 'nope.mdlus'"):
        manager.mark_best(job_id, "nope.mdlus")
    with pytest.raises(ValueError, match="no checkpoint"):
        manager.resolve_checkpoint(job_id, "../etc/passwd")


def test_bad_ids_and_specs_are_named_errors(manager, tmp_path):
    with pytest.raises(ValueError, match="invalid job id"):
        manager.status("../x")
    with pytest.raises(ValueError, match="no job 'nope'"):
        manager.status("nope")
    with pytest.raises(ValueError, match="unknown key 'Bogus'"):
        manager.start({**_spec(tmp_path), "Bogus": 1})
    assert not os.path.isdir(manager.runs_dir) or not os.listdir(manager.runs_dir)
    assert manager.list_jobs() == []
