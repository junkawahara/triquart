#pragma once

// Figure 7 of soda13.pdf with the Section 5 compression operations inserted at
// their specified points. Explicit frames avoid a recursive C++ call stack on
// caterpillars. Light subtrees get independent extracted component trees; the
// heavy child continues on the current HDT, contracting only after a size drop.
#include "traversal_stats.hpp"

#include <memory>
#include <utility>

namespace soda2013::detail {

struct LightTask {
    NodeId child;
    ComponentTree extracted;
};

template<class Decoration>
struct Frame {
    NodeId vertex;
    std::unique_ptr<Decoration> hdt;
    bool entered = false;
    std::vector<LightTask> light;
    std::size_t next_light = 0;
};

template<class Decoration>
class CountTraversal {
public:
    CountTraversal(const Tree& first, const Tree& second, CountTraversalStats& stats)
        : first_(first), second_(second), stats_(stats), matching_(first.size(), no_node),
          atom_for_original_(second.size(), no_node) {
        for (NodeId leaf : first.leaves()) matching_[leaf] = second.leaf(first.node(leaf).label);
    }

    using Result = decltype(std::declval<Decoration>().query());

    Result run() {
        Result result;
        std::vector<Frame<Decoration>> frames;
        frames.push_back({first_.root(), build(ComponentTree::from_tree(second_)), false, {}, 0});
        while (!frames.empty()) {
            auto& frame = frames.back();
            if (frame.entered) {
                if (frame.next_light == frame.light.size()) {
                    frames.pop_back();
                } else {
                    auto& task = frame.light[frame.next_light++];
                    const NodeId child = task.child;
                    auto next_hdt = build(std::move(task.extracted));
                    frames.push_back({child, std::move(next_hdt), false, {}, 0});
                }
                continue;
            }
            const NodeId v = frame.vertex;
            frame.entered = true;
            const std::size_t active = static_cast<std::size_t>(first_.subtree_leaves(v));
            const std::size_t components = frame.hdt->components().size();
            stats_.maximum_components_per_active_leaf = std::max(
                stats_.maximum_components_per_active_leaf, (components + active - 1) / active);
            if (first_.is_leaf(v)) {
                frame.hdt->recolor_atoms({{atom_for_original_[matching_[v]], 0}});  // Figure 7, leaf case.
                ++stats_.recolored_leaves;
                retire(frame.hdt);
                continue;
            }

            const auto& children = first_.node(v).children;
            NodeId heavy = children.front();
            for (NodeId child : children)
                if (first_.subtree_leaves(child) > first_.subtree_leaves(heavy)) heavy = child;

            std::vector<Recoloring> changes;
            Color color = 2;
            for (NodeId child : children) {
                if (child == heavy) continue;
                for (NodeId atom : leaf_atoms(child)) changes.push_back({atom, color});
                ++color;
            }
            frame.hdt->recolor_atoms(changes);
            stats_.recolored_leaves += changes.size();
            const auto shared = frame.hdt->query();
            result.accumulate(shared);
            ++stats_.queries;

            // Section 5: extract just AFTER the query, with each ci still
            // colored i. Only its leaves remain incontractible in the copy.
            for (NodeId child : children) {
                if (child == heavy) continue;
                auto extracted = frame.hdt->extract_atoms(leaf_atoms(child), &stats_.compression);
                frame.light.push_back({child, std::move(extracted)});
            }
            for (auto& change : changes) change.color = 0;
            frame.hdt->recolor_atoms(changes);
            stats_.recolored_leaves += changes.size();

            // A terminal greedy contraction has at most 8k atoms (see
            // component_tree.cpp). At this threshold each rebuild halves the
            // current size. Otherwise size < 16(k+1) <= 32k at Count(c1).
            const std::size_t k = static_cast<std::size_t>(first_.subtree_leaves(heavy));
            if (k < frame.hdt->components().size() / 16) {
                auto contracted = frame.hdt->contract_atoms(leaf_atoms(heavy), &stats_.compression);
                retire(frame.hdt);
                frame.hdt = build(std::move(contracted));
                ++stats_.heavy_contractions;
            }
            auto heavy_hdt = std::move(frame.hdt);
            frames.push_back({heavy, std::move(heavy_hdt), false, {}, 0});
        }
        return result;
    }

private:
    const Tree& first_;
    const Tree& second_;
    CountTraversalStats& stats_;
    std::vector<NodeId> matching_;
    // Section 5 paired pointers: only ONE decorated HDT is live in this DFS.
    // Pending light tasks hold component trees, so rebuilding these pointers
    // takes O(current size) without invalidating a suspended decoration.
    std::vector<NodeId> atom_for_original_;

    std::vector<NodeId> leaf_atoms(NodeId v) const {
        std::vector<NodeId> result;
        result.reserve(static_cast<std::size_t>(first_.subtree_leaves(v)));
        for (std::size_t i = first_.leaf_begin(v); i < first_.leaf_end(v); ++i)
            result.push_back(atom_for_original_[matching_[first_.leaves()[i]]]);
        return result;
    }
    std::unique_ptr<Decoration> build(ComponentTree tree) {
        ++stats_.hdt_builds;
        stats_.hdt_input_components += tree.size();
        auto result = std::make_unique<Decoration>(std::move(tree));
        for (NodeId v = 0; v < result->components().size(); ++v) {
            const NodeId original = result->components().atom(v).active_leaf;
            if (original != no_node) {
                atom_for_original_[original] = v;
                ++stats_.paired_leaf_updates;
            }
        }
        stats_.hdt_round_visits += result->hdt().stats().visited_components;
        return result;
    }
    void retire(std::unique_ptr<Decoration>& hdt) {
        const auto& updates = hdt->stats();
        stats_.updates.batches += updates.batches;
        stats_.updates.recomputed_components += updates.recomputed_components;
        stats_.updates.color_rows += updates.color_rows;
        stats_.updates.color_pairs += updates.color_pairs;
        hdt.reset();
    }
};

}  // namespace soda2013::detail
