#pragma once

#include "tree.hpp"

#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

#define CHECK(expression) do { if (!(expression)) \
    throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                             ": " #expression); } while (false)

namespace test_support {

template<class Exception, class Function>
void expect_exception(Function function) {
    bool caught = false;
    try { function(); } catch (const Exception&) { caught = true; }
    CHECK(caught);
}

inline soda2013::Tree random_tree(std::size_t leaves, std::mt19937_64& rng,
                                 std::size_t max_degree = 12) {
    using namespace soda2013;
    std::vector<Node> nodes;
    std::vector<NodeId> roots;
    for (std::size_t i = 0; i < leaves; ++i) {
        roots.push_back(nodes.size());
        nodes.push_back({std::to_string(i), {}});
    }
    while (roots.size() > 1) {
        std::shuffle(roots.begin(), roots.end(), rng);
        const auto degree = std::uniform_int_distribution<std::size_t>(
            2, std::min(max_degree, roots.size()))(rng);
        std::vector<NodeId> children(roots.end() - static_cast<std::ptrdiff_t>(degree), roots.end());
        roots.resize(roots.size() - degree);
        roots.push_back(nodes.size());
        nodes.push_back({{}, std::move(children)});
    }
    return Tree(std::move(nodes), roots.at(0));
}

inline soda2013::Tree caterpillar(std::size_t leaves) {
    using namespace soda2013;
    std::vector<Node> nodes;
    nodes.push_back({"0", {}});
    NodeId root = 0;
    for (std::size_t i = 1; i < leaves; ++i) {
        const NodeId leaf = nodes.size();
        nodes.push_back({std::to_string(i), {}});
        nodes.push_back({{}, {root, leaf}});
        root = nodes.size() - 1;
    }
    return Tree(std::move(nodes), root);
}

inline soda2013::Tree star(std::size_t leaves) {
    using namespace soda2013;
    if (leaves == 1) return Tree({{"0", {}}}, 0);
    std::vector<Node> nodes(1);
    for (std::size_t i = 0; i < leaves; ++i) {
        nodes[0].children.push_back(nodes.size());
        nodes.push_back({std::to_string(i), {}});
    }
    return Tree(std::move(nodes), 0);
}

inline soda2013::Tree balanced(std::size_t leaves) {
    using namespace soda2013;
    std::vector<Node> nodes;
    std::vector<NodeId> level;
    for (std::size_t i = 0; i < leaves; ++i) {
        level.push_back(nodes.size());
        nodes.push_back({std::to_string(i), {}});
    }
    while (level.size() > 1) {
        std::vector<NodeId> next;
        for (std::size_t i = 0; i < level.size(); i += 2) {
            if (i + 1 == level.size()) next.push_back(level[i]);
            else {
                next.push_back(nodes.size());
                nodes.push_back({{}, {level[i], level[i + 1]}});
            }
        }
        level = std::move(next);
    }
    return Tree(std::move(nodes), level.at(0));
}

}  // namespace test_support
