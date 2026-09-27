//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file operations/periodic.hpp
 * @brief Match the nodes of two boundary regions that a transform maps onto
 * each other: the node pairs a periodic boundary condition ties together.
 *
 * `match_periodic_nodes(mesh, slave, master, {transform})` maps every node of
 * the *slave* region through the affine transform and finds the *master*
 * node it lands on, within `mAtol`. A translation pairs opposite faces of a
 * box; a rotation pairs the cut faces of a sector. The result is the list a
 * Kratos periodic condition, or a Gmsh `$Periodic` section, is written from.
 *
 * ### Which nodes
 *
 * A region contributes the nodes it names: a Point region its entries, a Cell
 * region the nodes of its cells, a Side region the nodes of its facets. A node
 * of the slave region that the transform maps (within `mAtol`) onto itself,
 * and that is also in the master region -- a node on a rotation axis, say --
 * is a *fixed point*: it is counted, not paired.
 *
 * ### Guarantees
 *
 * The nearest master node within `mAtol` wins, ties to the lower node id. Two
 * slave nodes landing on one master node is an error (the tolerance is too
 * loose, or the regions overlap), and so is an unmatched slave node when
 * `mRequireComplete` is set. The pairs are ordered by slave node id, and the
 * result is the same on every backend and thread count.
 */

// System includes
#include <cstdint>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/region_ops.hpp"
#include "meshioplusplus/operations/transform.hpp"

namespace meshioplusplus {

/// How `match_periodic_nodes` matches.
struct PeriodicOptions {
    /// Maps a slave node's position onto its master's.
    AffineTransform mTransform;
    /// The largest distance, after the transform, at which two nodes match.
    /// Must be positive.
    double mAtol = 1e-8;
    /// Throw when a slave node has no master within `mAtol`.
    bool mRequireComplete = true;
};

/// The matched node pairs.
struct PeriodicPairs {
    /// Int64 `(k,)`: the slave node ids, ascending.
    NDArray mSlave;
    /// Int64 `(k,)`: the master node each slave node maps onto.
    NDArray mMaster;
    /// Int64: slave nodes with no master within `mAtol` (empty when
    /// `mRequireComplete` is set, since those throw).
    NDArray mUnmatched;
    /// Slave nodes the transform leaves in place that are also master nodes.
    std::int64_t mNumFixed = 0;
    /// The largest distance between a transformed slave node and its master.
    double mMaxResidual = 0.0;
};

/**
 * @brief The master node each node of the @p rSlave region maps onto.
 * @throws std::invalid_argument on a missing or ambiguous region, a
 *         non-positive tolerance, a master node claimed twice, or (with
 *         `mRequireComplete`) an unmatched slave node.
 */
MESHIOPLUSPLUS_API PeriodicPairs match_periodic_nodes(const Mesh& rMesh,
                                                      const RegionSelector& rSlave,
                                                      const RegionSelector& rMaster,
                                                      const PeriodicOptions& rOptions = {});

}  // namespace meshioplusplus
