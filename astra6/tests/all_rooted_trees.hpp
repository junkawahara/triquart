#pragma once

// Test-only exhaustive set partitions enumerate all labelled rooted trees.
#include <functional>
#include <string>
#include <vector>

inline std::vector<std::string> all_subtrees(const std::vector<std::size_t>& labels) {
    if (labels.size() == 1) return {std::to_string(labels[0])};
    std::vector<std::string> result;
    std::vector<std::vector<std::size_t>> blocks;
    std::function<void(std::size_t)> partition = [&](std::size_t pos) {
        if (pos == labels.size()) {
            if (blocks.size() < 2) return;
            std::vector<std::string> combinations{"("};
            for (std::size_t i = 0; i < blocks.size(); ++i) {
                const auto subtrees = all_subtrees(blocks[i]);
                std::vector<std::string> next;
                for (const auto& prefix : combinations)
                    for (const auto& subtree : subtrees) next.push_back(prefix + (i ? "," : "") + subtree);
                combinations = std::move(next);
            }
            for (auto& tree : combinations) result.push_back(tree + ")");
            return;
        }
        // Do not retain references to blocks elements over recursive pushes.
        const std::size_t count = blocks.size();
        for (std::size_t i = 0; i < count; ++i) {
            blocks[i].push_back(labels[pos]);
            partition(pos + 1);
            blocks[i].pop_back();
        }
        blocks.push_back({labels[pos]});
        partition(pos + 1);
        blocks.pop_back();
    };
    partition(0);
    return result;
}
