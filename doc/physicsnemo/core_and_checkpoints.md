---
title: Core and checkpoints
description: physicsnemo.Module, the .mdlus format, and the model card that gives its tensors meaning.
---

# Core and checkpoints

`physicsnemo.core` is small and you will use two things from it constantly: `Module` and the checkpoint format it defines.

## `physicsnemo.Module`

A `physicsnemo.Module` **is** a `torch.nn.Module` — it subclasses it, so everything you know about `forward`, `parameters()`, `.to(device)` and `state_dict()` applies unchanged. What it adds is that it **records the arguments it was constructed with**.

That single addition is the reason the checkpoint format works. A plain torch `state_dict` is a bag of tensors: to load it you must first reconstruct the architecture yourself, in code, with exactly the hyperparameters it was trained with. A `physicsnemo.Module` writes those hyperparameters into the checkpoint, so loading needs no architecture code at all:

```python
import physicsnemo
model = physicsnemo.Module.from_checkpoint("surrogate.mdlus")   # that is all
```

`ModelMetaData` is the descriptor a model class carries (name, whether it supports CUDA graphs, AMP, ONNX export, and so on). You only touch it when writing a new architecture.

## What is inside a `.mdlus` file

A zip archive containing:

- the `state_dict` — the weights;
- the constructor arguments — the architecture;
- the class's registry name, so `from_checkpoint` can find the class again.

![A .mdlus archive, a TorchScript file and a model-card sidecar side by side](/images/physicsnemo/mdlus_and_card.svg)

*Figure 1: The two checkpoint formats, the card that gives their tensors meaning, and what a deployment does with all three every step.*

The consequences matter in practice:

- **A `.mdlus` needs physicsnemo installed to load.** The class it names has to exist. If you want an artifact that runs without physicsnemo, export to ONNX (see [Diffusion and deployment](./diffusion_and_deployment.md)).
- **A `.mdlus` does not need your training script.** This is the point.
- **It does not record what the numbers mean.** Nothing in the format says which field is in channel 0, what units it is in, or whether the targets were normalized. That gap is what the model card below fills.

`physicsnemo.core.registry` is the class registry `from_checkpoint` looks names up in; `physicsnemo.compat` handles loading checkpoints written by older releases.

## TorchScript, the other format

A TorchScript file (`.pt` written by `torch.jit.script`/`trace`) is the alternative: it stores the *traced computation*, not the architecture, so it loads with only torch installed and no physicsnemo at all. In exchange it is frozen — no `torch.compile`, no easy surgery, and anything not traceable (a gpytorch head, for instance) simply cannot go in one.

Rule of thumb: **`.mdlus` while you are still training, TorchScript or ONNX when you ship.**

## The model card

A checkpoint that says nothing about its fields is a checkpoint you can silently misuse. Writing a model's raw output onto a temperature field when it was trained on targets scaled to zero mean and unit variance produces finite, plausible-looking, completely wrong numbers.

So the useful convention is a **model card**: a JSON sidecar next to the checkpoint naming the input and output fields, their order, and — the part that prevents the failure above — the scaling the targets were trained under, plus its mirror for the features. Anything that loads the checkpoint validates its configuration against the card, standardizes what it feeds the model, and inverts the output normalization before touching a field. An absent key is exactly the identity, so cards stay optional; a wrong channel count raises rather than producing silent NaNs.

This is the single most useful thing to know about deploying a surrogate: **if you train a model, save it with a card.**

## In meshio++

`meshioplusplus.physicsnemo.run_training` writes both. Every run directory holds `best.mdlus` and `last.mdlus` written through `physicsnemo.Module.save`, plus a `*.card.json` recording the feature and target column names in order, the per-column normalization actually applied, the graph or grid construction rule, and the schema version those columns were produced under.

`predict` reads the card back, checks the mesh it is handed still produces the same columns, and raises by name when it does not — a mesh missing a field the model was trained on is a named error, never a quietly shorter feature matrix. See [PhysicsNeMo integration](../physicsnemo.md#training-and-prediction).

Next: [Models](./models.md) — what to put in the checkpoint.
