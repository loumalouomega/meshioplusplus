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
// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

// Project includes
#include "weighted_matching.hpp"

namespace meshioplusplus::detail {

namespace {

using I = std::int64_t;

// The state of one run. Vertices are 0..n-1, blossoms n..2n-1. An edge k
// has the endpoints 2k (its first vertex) and 2k + 1 (its second); `p ^ 1`
// is the other end of endpoint p. Labels: 0 free, 1 S (outer), 2 T
// (inner); 4 marks a blossom on the path being scanned (5 = S marked).
class WmState {
public:
    WmState(I N, const std::vector<MatchingEdge>& rEdges, bool MaxCardinality)
        : n(N), mMaxCardinality(MaxCardinality) {
        I most = 0;
        for (const MatchingEdge& r_e : rEdges) {
            if (r_e.mA == r_e.mB)
                continue;
            if (r_e.mA < 0 || r_e.mB < 0 || r_e.mA >= N || r_e.mB >= N)
                throw std::invalid_argument("meshio++: max_weight_matching: vertex out of range");
            // Doubled: every slack and dual stays an integer.
            mEdges.push_back({r_e.mA, r_e.mB, 2 * r_e.mWeight});
            most = std::max(most, 2 * r_e.mWeight);
        }
        const auto m = static_cast<I>(mEdges.size());
        mEndpoint.resize(static_cast<std::size_t>(2 * m));
        mNeighbend.assign(static_cast<std::size_t>(n), {});
        for (I k = 0; k < m; ++k) {
            const MatchingEdge& r_e = mEdges[static_cast<std::size_t>(k)];
            mEndpoint[static_cast<std::size_t>(2 * k)] = r_e.mA;
            mEndpoint[static_cast<std::size_t>(2 * k + 1)] = r_e.mB;
            mNeighbend[static_cast<std::size_t>(r_e.mA)].push_back(2 * k + 1);
            mNeighbend[static_cast<std::size_t>(r_e.mB)].push_back(2 * k);
        }
        const auto nn = static_cast<std::size_t>(2 * n);
        mMate.assign(static_cast<std::size_t>(n), -1);
        mLabel.assign(nn, 0);
        mLabelEnd.assign(nn, -1);
        mInBlossom.resize(static_cast<std::size_t>(n));
        for (I v = 0; v < n; ++v)
            mInBlossom[static_cast<std::size_t>(v)] = v;
        mBlossomParent.assign(nn, -1);
        mBlossomChilds.assign(nn, {});
        mBlossomBase.assign(nn, -1);
        for (I v = 0; v < n; ++v)
            mBlossomBase[static_cast<std::size_t>(v)] = v;
        mBlossomEndps.assign(nn, {});
        mBestEdge.assign(nn, -1);
        mBlossomBestEdges.assign(nn, {});
        mHasBestEdges.assign(nn, 0);
        for (I b = 2 * n - 1; b >= n; --b)
            mUnused.push_back(b);
        mDual.assign(nn, 0);
        for (I v = 0; v < n; ++v)
            mDual[static_cast<std::size_t>(v)] = most;
        mAllowEdge.assign(static_cast<std::size_t>(m), 0);
    }

    std::vector<I> Run() {
        for (I stage = 0; stage < n; ++stage) {
            std::fill(mLabel.begin(), mLabel.end(), 0);
            std::fill(mBestEdge.begin(), mBestEdge.end(), -1);
            for (I b = n; b < 2 * n; ++b) {
                mBlossomBestEdges[static_cast<std::size_t>(b)].clear();
                mHasBestEdges[static_cast<std::size_t>(b)] = 0;
            }
            std::fill(mAllowEdge.begin(), mAllowEdge.end(), 0);
            mQueue.clear();
            for (I v = 0; v < n; ++v)
                if (Mate(v) == -1 && Label(InBlossom(v)) == 0)
                    AssignLabel(v, 1, -1);
            bool augmented = false;
            for (;;) {
                while (!mQueue.empty() && !augmented) {
                    const I v = mQueue.back();
                    mQueue.pop_back();
                    for (const I p : mNeighbend[static_cast<std::size_t>(v)]) {
                        const I k = p / 2;
                        const I w = Endpoint(p);
                        if (InBlossom(v) == InBlossom(w))
                            continue;
                        I kslack = 0;
                        if (!mAllowEdge[static_cast<std::size_t>(k)]) {
                            kslack = Slack(k);
                            if (kslack <= 0)
                                mAllowEdge[static_cast<std::size_t>(k)] = 1;
                        }
                        if (mAllowEdge[static_cast<std::size_t>(k)]) {
                            if (Label(InBlossom(w)) == 0) {
                                AssignLabel(w, 2, p ^ 1);
                            } else if (Label(InBlossom(w)) == 1) {
                                const I base = ScanBlossom(v, w);
                                if (base >= 0) {
                                    AddBlossom(base, k);
                                } else {
                                    AugmentMatching(k);
                                    augmented = true;
                                    break;
                                }
                            } else if (Label(w) == 0) {
                                Label(w) = 2;
                                LabelEnd(w) = p ^ 1;
                            }
                        } else if (Label(InBlossom(w)) == 1) {
                            const I b = InBlossom(v);
                            if (BestEdge(b) == -1 || kslack < Slack(BestEdge(b)))
                                BestEdge(b) = k;
                        } else if (Label(w) == 0) {
                            if (BestEdge(w) == -1 || kslack < Slack(BestEdge(w)))
                                BestEdge(w) = k;
                        }
                    }
                }
                if (augmented)
                    break;
                // The least of the four slacks.
                int type = -1;
                I delta = 0, delta_edge = -1, delta_blossom = -1;
                if (!mMaxCardinality) {
                    type = 1;
                    delta = *std::min_element(mDual.begin(), mDual.begin() + n);
                }
                for (I v = 0; v < n; ++v)
                    if (Label(InBlossom(v)) == 0 && BestEdge(v) != -1) {
                        const I d = Slack(BestEdge(v));
                        if (type == -1 || d < delta) {
                            delta = d;
                            type = 2;
                            delta_edge = BestEdge(v);
                        }
                    }
                for (I b = 0; b < 2 * n; ++b)
                    if (Parent(b) == -1 && Label(b) == 1 && BestEdge(b) != -1) {
                        const I d = Slack(BestEdge(b)) / 2;
                        if (type == -1 || d < delta) {
                            delta = d;
                            type = 3;
                            delta_edge = BestEdge(b);
                        }
                    }
                for (I b = n; b < 2 * n; ++b)
                    if (Base(b) >= 0 && Parent(b) == -1 && Label(b) == 2 &&
                        (type == -1 || Dual(b) < delta)) {
                        delta = Dual(b);
                        type = 4;
                        delta_blossom = b;
                    }
                if (type == -1) {
                    // Max cardinality and no further progress: stop the run.
                    type = 1;
                    delta = std::max<I>(0, *std::min_element(mDual.begin(), mDual.begin() + n));
                }
                for (I v = 0; v < n; ++v) {
                    const I l = Label(InBlossom(v));
                    if (l == 1)
                        Dual(v) -= delta;
                    else if (l == 2)
                        Dual(v) += delta;
                }
                for (I b = n; b < 2 * n; ++b)
                    if (Base(b) >= 0 && Parent(b) == -1) {
                        if (Label(b) == 1)
                            Dual(b) += delta;
                        else if (Label(b) == 2)
                            Dual(b) -= delta;
                    }
                if (type == 1)
                    break;
                if (type == 2) {
                    mAllowEdge[static_cast<std::size_t>(delta_edge)] = 1;
                    I i = Edge(delta_edge).mA, j = Edge(delta_edge).mB;
                    if (Label(InBlossom(i)) == 0)
                        std::swap(i, j);
                    mQueue.push_back(i);
                } else if (type == 3) {
                    mAllowEdge[static_cast<std::size_t>(delta_edge)] = 1;
                    mQueue.push_back(Edge(delta_edge).mA);
                } else {
                    ExpandBlossom(delta_blossom, false);
                }
            }
            if (!augmented)
                break;
            for (I b = n; b < 2 * n; ++b)
                if (Parent(b) == -1 && Base(b) >= 0 && Label(b) == 1 && Dual(b) == 0)
                    ExpandBlossom(b, true);
        }
        std::vector<I> out(static_cast<std::size_t>(n), -1);
        for (I v = 0; v < n; ++v)
            if (Mate(v) >= 0)
                out[static_cast<std::size_t>(v)] = Endpoint(Mate(v));
        return out;
    }

private:
    I n = 0;
    bool mMaxCardinality = false;
    std::vector<MatchingEdge> mEdges;
    std::vector<I> mEndpoint;
    std::vector<std::vector<I>> mNeighbend;
    std::vector<I> mMate, mLabel, mLabelEnd, mInBlossom, mBlossomParent, mBlossomBase;
    std::vector<std::vector<I>> mBlossomChilds, mBlossomEndps, mBlossomBestEdges;
    std::vector<char> mHasBestEdges;
    std::vector<I> mBestEdge, mUnused, mDual, mQueue;
    std::vector<char> mAllowEdge;

    static std::size_t Z(I X) { return static_cast<std::size_t>(X); }
    I& Mate(I V) { return mMate[Z(V)]; }
    I& Label(I B) { return mLabel[Z(B)]; }
    I& LabelEnd(I B) { return mLabelEnd[Z(B)]; }
    I& InBlossom(I V) { return mInBlossom[Z(V)]; }
    I& Parent(I B) { return mBlossomParent[Z(B)]; }
    I& Base(I B) { return mBlossomBase[Z(B)]; }
    I& BestEdge(I B) { return mBestEdge[Z(B)]; }
    I& Dual(I B) { return mDual[Z(B)]; }
    I Endpoint(I P) const { return mEndpoint[Z(P)]; }
    const MatchingEdge& Edge(I K) const { return mEdges[Z(K)]; }

    I Slack(I K) const {
        const MatchingEdge& r_e = mEdges[Z(K)];
        return mDual[Z(r_e.mA)] + mDual[Z(r_e.mB)] - 2 * r_e.mWeight;
    }

    void Leaves(I B, std::vector<I>& rOut) const {
        if (B < n) {
            rOut.push_back(B);
            return;
        }
        for (const I c : mBlossomChilds[Z(B)])
            Leaves(c, rOut);
    }

    std::vector<I> Leaves(I B) const {
        std::vector<I> out;
        Leaves(B, out);
        return out;
    }

    // Label the top blossom of W with T, and the blossom across its base's
    // mate with S, as the tree grows through P.
    void AssignLabel(I W, I T, I P) {
        const I b = InBlossom(W);
        Label(W) = Label(b) = T;
        LabelEnd(W) = LabelEnd(b) = P;
        BestEdge(W) = BestEdge(b) = -1;
        if (T == 1) {
            Leaves(b, mQueue);
        } else if (T == 2) {
            const I base = Base(b);
            AssignLabel(Endpoint(Mate(base)), 1, Mate(base) ^ 1);
        }
    }

    // Trace back from V and W to the trees' roots: a new blossom's base, or
    // -1 when the roots differ (an augmenting path).
    I ScanBlossom(I V, I W) {
        std::vector<I> path;
        I base = -1;
        while (V != -1 || W != -1) {
            I b = InBlossom(V);
            if (Label(b) & 4) {
                base = Base(b);
                break;
            }
            path.push_back(b);
            Label(b) = 5;
            if (LabelEnd(b) == -1) {
                V = -1;
            } else {
                V = Endpoint(LabelEnd(b));
                b = InBlossom(V);
                V = Endpoint(LabelEnd(b));
            }
            if (W != -1)
                std::swap(V, W);
        }
        for (const I b : path)
            Label(b) = 1;
        return base;
    }

    // Shrink the cycle closed by edge K through the trees from base Base.
    void AddBlossom(I BaseV, I K) {
        I v = Edge(K).mA, w = Edge(K).mB;
        const I bb = InBlossom(BaseV);
        I bv = InBlossom(v), bw = InBlossom(w);
        const I b = mUnused.back();
        mUnused.pop_back();
        Base(b) = BaseV;
        Parent(b) = -1;
        Parent(bb) = b;
        std::vector<I>& r_path = mBlossomChilds[Z(b)];
        std::vector<I>& r_endps = mBlossomEndps[Z(b)];
        r_path.clear();
        r_endps.clear();
        while (bv != bb) {
            Parent(bv) = b;
            r_path.push_back(bv);
            r_endps.push_back(LabelEnd(bv));
            v = Endpoint(LabelEnd(bv));
            bv = InBlossom(v);
        }
        r_path.push_back(bb);
        std::reverse(r_path.begin(), r_path.end());
        std::reverse(r_endps.begin(), r_endps.end());
        r_endps.push_back(2 * K);
        while (bw != bb) {
            Parent(bw) = b;
            r_path.push_back(bw);
            r_endps.push_back(LabelEnd(bw) ^ 1);
            w = Endpoint(LabelEnd(bw));
            bw = InBlossom(w);
        }
        Label(b) = 1;
        LabelEnd(b) = LabelEnd(bb);
        Dual(b) = 0;
        for (const I leaf : Leaves(b)) {
            if (Label(InBlossom(leaf)) == 2)
                mQueue.push_back(leaf);
            InBlossom(leaf) = b;
        }
        // The least-slack edge from the new blossom to each S-blossom.
        std::vector<I> best_to(Z(2 * n), -1);
        for (const I sub : mBlossomChilds[Z(b)]) {
            std::vector<std::vector<I>> lists;
            if (!mHasBestEdges[Z(sub)]) {
                for (const I leaf : Leaves(sub)) {
                    std::vector<I> l;
                    for (const I p : mNeighbend[Z(leaf)])
                        l.push_back(p / 2);
                    lists.push_back(std::move(l));
                }
            } else {
                lists.push_back(mBlossomBestEdges[Z(sub)]);
            }
            for (const auto& r_list : lists)
                for (const I k : r_list) {
                    I i = Edge(k).mA, j = Edge(k).mB;
                    if (InBlossom(j) == b)
                        std::swap(i, j);
                    const I bj = InBlossom(j);
                    if (bj != b && Label(bj) == 1 &&
                        (best_to[Z(bj)] == -1 || Slack(k) < Slack(best_to[Z(bj)])))
                        best_to[Z(bj)] = k;
                }
            mBlossomBestEdges[Z(sub)].clear();
            mHasBestEdges[Z(sub)] = 0;
            BestEdge(sub) = -1;
        }
        std::vector<I>& r_best = mBlossomBestEdges[Z(b)];
        r_best.clear();
        for (const I k : best_to)
            if (k != -1)
                r_best.push_back(k);
        mHasBestEdges[Z(b)] = 1;
        BestEdge(b) = -1;
        for (const I k : r_best)
            if (BestEdge(b) == -1 || Slack(k) < Slack(BestEdge(b)))
                BestEdge(b) = k;
    }

    // A child's index from the end when negative.
    static std::size_t Wrap(I J, std::size_t N) {
        const I m = static_cast<I>(N);
        return static_cast<std::size_t>(((J % m) + m) % m);
    }

    // Expand blossom B into its children; at a stage's end recursively
    // those whose dual is zero, otherwise relabel the T-blossom's children.
    void ExpandBlossom(I B, bool EndStage) {
        const std::vector<I> childs = mBlossomChilds[Z(B)];
        for (const I s : childs) {
            Parent(s) = -1;
            if (s < n) {
                InBlossom(s) = s;
            } else if (EndStage && Dual(s) == 0) {
                ExpandBlossom(s, EndStage);
            } else {
                for (const I leaf : Leaves(s))
                    InBlossom(leaf) = s;
            }
        }
        if (!EndStage && Label(B) == 2) {
            const std::vector<I>& r_endps = mBlossomEndps[Z(B)];
            const std::size_t nc = childs.size();
            const I entry = InBlossom(Endpoint(LabelEnd(B) ^ 1));
            I j = static_cast<I>(std::find(childs.begin(), childs.end(), entry) - childs.begin());
            I step = 0, trick = 0;
            if (j & 1) {
                j -= static_cast<I>(nc);
                step = 1;
                trick = 0;
            } else {
                step = -1;
                trick = 1;
            }
            I p = LabelEnd(B);
            while (j != 0) {
                Label(Endpoint(p ^ 1)) = 0;
                Label(Endpoint(r_endps[Wrap(j - trick, nc)] ^ trick ^ 1)) = 0;
                AssignLabel(Endpoint(p ^ 1), 2, p);
                mAllowEdge[Z(r_endps[Wrap(j - trick, nc)] / 2)] = 1;
                j += step;
                p = r_endps[Wrap(j - trick, nc)] ^ trick;
                mAllowEdge[Z(p / 2)] = 1;
                j += step;
            }
            const I bv = childs[Wrap(j, nc)];
            Label(Endpoint(p ^ 1)) = Label(bv) = 2;
            LabelEnd(Endpoint(p ^ 1)) = LabelEnd(bv) = p;
            BestEdge(bv) = -1;
            j += step;
            while (childs[Wrap(j, nc)] != entry) {
                const I sub = childs[Wrap(j, nc)];
                if (Label(sub) == 1) {
                    j += step;
                    continue;
                }
                I found = -1;
                for (const I leaf : Leaves(sub))
                    if (Label(leaf) != 0) {
                        found = leaf;
                        break;
                    }
                if (found >= 0) {
                    Label(found) = 0;
                    Label(Endpoint(Mate(Base(sub)))) = 0;
                    AssignLabel(found, 2, LabelEnd(found));
                }
                j += step;
            }
        }
        Label(B) = LabelEnd(B) = -1;
        mBlossomChilds[Z(B)].clear();
        mBlossomEndps[Z(B)].clear();
        Base(B) = -1;
        mBlossomBestEdges[Z(B)].clear();
        mHasBestEdges[Z(B)] = 0;
        BestEdge(B) = -1;
        mUnused.push_back(B);
    }

    // Swap matched and unmatched edges along the even path from V to B's
    // base, rotating B's children so V's child is its new base.
    void AugmentBlossom(I B, I V) {
        I t = V;
        while (Parent(t) != B)
            t = Parent(t);
        if (t >= n)
            AugmentBlossom(t, V);
        std::vector<I>& r_childs = mBlossomChilds[Z(B)];
        std::vector<I>& r_endps = mBlossomEndps[Z(B)];
        const std::size_t nc = r_childs.size();
        const I i =
            static_cast<I>(std::find(r_childs.begin(), r_childs.end(), t) - r_childs.begin());
        I j = i, step = 0, trick = 0;
        if (i & 1) {
            j -= static_cast<I>(nc);
            step = 1;
            trick = 0;
        } else {
            step = -1;
            trick = 1;
        }
        while (j != 0) {
            j += step;
            t = r_childs[Wrap(j, nc)];
            const I p = r_endps[Wrap(j - trick, nc)] ^ trick;
            if (t >= n)
                AugmentBlossom(t, Endpoint(p));
            j += step;
            t = r_childs[Wrap(j, nc)];
            if (t >= n)
                AugmentBlossom(t, Endpoint(p ^ 1));
            Mate(Endpoint(p)) = p ^ 1;
            Mate(Endpoint(p ^ 1)) = p;
        }
        std::rotate(r_childs.begin(), r_childs.begin() + i, r_childs.end());
        std::rotate(r_endps.begin(), r_endps.begin() + i, r_endps.end());
        Base(B) = Base(r_childs[0]);
    }

    // Augment along the path through edge K between two S-vertices.
    void AugmentMatching(I K) {
        const I v = Edge(K).mA, w = Edge(K).mB;
        for (int side = 0; side < 2; ++side) {
            I s = side == 0 ? v : w;
            I p = side == 0 ? 2 * K + 1 : 2 * K;
            for (;;) {
                const I bs = InBlossom(s);
                if (bs >= n)
                    AugmentBlossom(bs, s);
                Mate(s) = p;
                if (LabelEnd(bs) == -1)
                    break;
                const I t = Endpoint(LabelEnd(bs));
                const I bt = InBlossom(t);
                s = Endpoint(LabelEnd(bt));
                const I j = Endpoint(LabelEnd(bt) ^ 1);
                if (bt >= n)
                    AugmentBlossom(bt, j);
                Mate(j) = LabelEnd(bt);
                p = LabelEnd(bt) ^ 1;
            }
        }
    }
};

}  // namespace

std::vector<std::int64_t> max_weight_matching(std::int64_t NumVertices,
                                              const std::vector<MatchingEdge>& rEdges,
                                              bool MaxCardinality) {
    if (NumVertices < 0)
        throw std::invalid_argument("meshio++: max_weight_matching: negative vertex count");
    if (NumVertices == 0)
        return {};
    WmState state(NumVertices, rEdges, MaxCardinality);
    return state.Run();
}

}  // namespace meshioplusplus::detail
