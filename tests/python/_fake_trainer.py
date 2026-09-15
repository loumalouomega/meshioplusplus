"""A stand-in for ``python -m meshioplusplus.physicsnemo.train`` used by the
job-manager tests (selected via ``MESHIOPLUSPLUS_TRAIN_COMMAND``): reads the
spec, writes the run-directory files the real trainer writes -- one metrics
row and a progress update per "epoch", a fake ``.mdlus`` plus card every
``CheckpointEvery`` epochs, ``final.mdlus`` at the end -- and honours SIGTERM
the way the real one does (finish the epoch, write final, exit 143).
``Notes == "fail"`` makes it exit 3 after two epochs, ``Notes == "hang"`` makes
it ignore SIGTERM (so the SIGKILL path is exercised).
"""

import argparse
import json
import os
import signal
import sys
import time

# NOTE: deliberately no sys.path manipulation. The manager spawns this with
# `sys.executable`, so the meshioplusplus the test session imported is the one
# importable here -- while putting `src/python` on the path would SHADOW an
# installed package with a source tree that has no compiled `_core` (which is
# exactly how CI installs it, so that shadowing would fail there and pass
# locally against an editable install).
from meshioplusplus.physicsnemo._train import (
    CHECKPOINT_DIR,
    PROGRESS_FILE,
    append_metrics,
    card_path,
    load_spec,
    metrics_row,
    write_json_atomic,
)

STOP = {"requested": False}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--epoch-seconds", type=float, default=0.05)
    args = ap.parse_args()
    spec = load_spec(args.spec)
    run_dir = spec.resolved_run_dir()
    if spec.notes != "hang":
        signal.signal(signal.SIGTERM, lambda *_: STOP.update(requested=True))
    print(f"fake trainer: {spec.epochs} epochs on {spec.manifest}", flush=True)
    start = time.time()
    best = None
    for epoch in range(spec.epochs):
        time.sleep(args.epoch_seconds)
        train_loss = 1.0 / (epoch + 1)
        valid_loss = 1.2 / (epoch + 1)
        append_metrics(
            run_dir,
            metrics_row(
                epoch,
                train_loss,
                valid_loss,
                spec.learning_rate,
                time.time() - start,
                args.epoch_seconds,
            ),
        )
        if best is None or valid_loss < best[1]:
            best = (epoch, valid_loss)
        write_json_atomic(
            os.path.join(run_dir, PROGRESS_FILE),
            {
                "epoch": epoch + 1,
                "epochs": spec.epochs,
                "best_epoch": best[0],
                "best_valid_loss": best[1],
                "best_checkpoint": None,
                "eta_seconds": (spec.epochs - epoch - 1) * args.epoch_seconds,
                "device": "cpu",
                "completed": False,
            },
        )
        if spec.checkpoint_every and (epoch + 1) % spec.checkpoint_every == 0:
            _checkpoint(run_dir, f"Fake.0.{epoch}.mdlus", epoch, valid_loss)
        print(
            f"epoch {epoch} train {train_loss:.3e} valid {valid_loss:.3e}", flush=True
        )
        if spec.notes == "fail" and epoch >= 1:
            print("fake trainer: failing on purpose", flush=True)
            sys.exit(3)
        if STOP["requested"]:
            _checkpoint(run_dir, "final.mdlus", epoch, valid_loss)
            print("fake trainer: stopped", flush=True)
            sys.exit(143)
    _checkpoint(run_dir, "final.mdlus", spec.epochs - 1, best[1])
    _checkpoint(run_dir, "best.mdlus", best[0], best[1])
    progress = json.load(open(os.path.join(run_dir, PROGRESS_FILE)))
    progress["completed"] = True
    progress["best_checkpoint"] = os.path.join(run_dir, CHECKPOINT_DIR, "best.mdlus")
    write_json_atomic(os.path.join(run_dir, PROGRESS_FILE), progress)
    print("fake trainer: done", flush=True)


def _checkpoint(run_dir, name, epoch, valid_loss):
    path = os.path.join(run_dir, CHECKPOINT_DIR, name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(b"not a real checkpoint")
    write_json_atomic(
        card_path(path),
        {"version": 1, "epoch": epoch, "valid_loss": valid_loss, "checkpoint": name},
    )


if __name__ == "__main__":
    main()
