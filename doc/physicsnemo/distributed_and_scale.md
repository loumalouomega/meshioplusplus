---
title: Distributed and scale
description: DistributedManager, process groups, ShardTensor - two kinds of parallelism, and which of them one GPU can actually run.
---

# Distributed and scale

Two different subpackages, two different kinds of parallelism, and confusing them is the usual mistake.

![Data parallelism next to domain parallelism](/images/physicsnemo/parallelism.svg)

*Figure 1: Data parallelism replicates the model and splits the batch; domain parallelism splits the sample itself.*

## `physicsnemo.distributed` — data parallelism

Every rank holds the **whole model** and a **different slice of the data**. Gradients are averaged across ranks; the model is replicated.

- `DistributedManager` — the singleton that owns rank, world size, device and the backend. Everything else asks it.
- `ProcessGroupConfig`, `ProcessGroupNode` — subgroups of ranks, for models that want more than one axis of parallelism.
- Collectives with autograd support: `gather_v`, `all_gather_v`, `scatter_v`, `indexed_all_to_all_v`, `reduce_loss`, `fused_all_reduce`.
- `mark_module_as_shared` — tells the gradient machinery a module's parameters are replicated rather than sharded.
- `distributed.fft` — distributed FFTs, for spectral models.

### The alignment problem

A solver has its own communicator abstraction, and torch has `torch.distributed`. If the two disagree about which rank you are, the failure is silent and the results are garbage.

The fix is not subtle — align them at start-up and **check loudly** when they disagree — but it has to be done deliberately, because nothing in either library will notice on its own.

### Halo partitioning

Splitting a mesh graph across ranks naively truncates the interfaces: a node near a partition boundary loses the neighbours that live on another rank, so its one-hop neighbourhood no longer matches a serial run and the message-passing result is wrong at exactly the places that matter.

The remedy is per-rank subgraphs whose owned sets partition the global node set exactly, and whose one-hop neighbourhoods do match — a halo of boundary cells around each partition, with gradients synced by `DistributedDataParallel` on top.

## `physicsnemo.domain_parallel` — domain parallelism

Every rank holds **part of one tensor**. The mesh itself is split, not the batch. This is what you want when a single sample does not fit on one GPU.

- `ShardTensor` — a tensor sharded across a device mesh, with `shard_tensor`, `scatter_tensor` and `sync_module_over_mesh`.

Note it moved: `ShardTensor` lives in `physicsnemo.domain_parallel`, **not** in `physicsnemo.distributed`, as of 2.2.

Two facts about what a CPU-only or single-GPU machine can test here, because they are easy to get wrong in both directions. `ShardTensor` itself needs only a torch `DeviceMesh`, which `init_device_mesh("cpu", ...)` provides, and upstream's operator handlers register on a CUDA-less host once asked to — so the *semantics* are testable over gloo. What is not testable is the NCCL transport, plus a handful of genuinely CUDA-only remnants: `sharding_shapes="infer"`, ring attention, and the sharded kNN and radius search.

The trap worth pinning: **the backward of a mean over a `ShardTensor` already leaves the serial gradient on every rank**, so a DDP-style all-reduce on top is wrong by the rank count. It produces a plausible number, which is why it needs an assertion rather than a reading.

## Sharded checkpoints

An FSDP2 model's parameters are `DTensor`s — each rank holds a shard. Saving them directly produces a checkpoint that reports success and cannot be loaded, with each rank having written only its own piece. Gather them to full tensors and write from rank 0.

## In meshio++

Nothing here is implemented, and the reason is worth stating precisely rather than as a gap: **meshio++ is a library, not a process.** It has no communicator, launches nothing, and has no opinion about ranks. What it has instead is the piece a distributed training run needs from a mesh library, which is the *decomposition itself*:

[`partition`](../partition.md) splits a mesh into exactly N balanced pieces — by a space-filling curve always, or by KaHIP's edge-cut minimization when it is built in — and `ghost_layers=N` grows each piece by that many shared-node breadth-first layers. That halo is exactly the one the graph-partitioning paragraph above describes, and shared-node rather than shared-facet adjacency is the conservative choice: it is what a node-based assembly needs, and it is a superset of what message passing on a mesh-edge graph needs.

`partition_labels` gives the same decomposition as a flat per-cell array instead of N meshes, which is the form to hand a distributed data loader. It **refuses** a nonzero `ghost_layers` rather than ignoring it, because a flat label array is the *ownership* map and a cell can be a ghost of several parts at once — there is no honest single label for it.

The rest — aligning two communicators, gathering DTensor shards, choosing NCCL over gloo — belongs to whatever process owns the training loop.

Next: [Companion packages](./companion_packages.md).
