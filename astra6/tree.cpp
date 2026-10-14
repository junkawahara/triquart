// Input representation and validation for soda13.pdf, Sections 2 and 4.
#include "tree.hpp"

#include <cctype>
#include <cmath>
#include <utility>

namespace soda2013 {

Tree::Tree(std::vector<Node> nodes, NodeId root)
    : nodes_(std::move(nodes)), root_(root), parents_(nodes_.size(), no_node),
      depths_(nodes_.size()), begin_(nodes_.size()), end_(nodes_.size()) {
    if (nodes_.empty() || root >= size())
        throw std::invalid_argument("tree needs a valid root");
    std::vector<bool> seen(size(), false);
    std::vector<NodeId> stack{root};
    seen[root] = true;
    while (!stack.empty()) {
        const NodeId v = stack.back();
        stack.pop_back();
        preorder_.push_back(v);
        const auto& children = nodes_[v].children;
        if (children.size() == 1)
            throw std::invalid_argument("unary internal nodes are not supported");
        if (children.empty()) {
            if (nodes_[v].label.empty() || !labels_.emplace(nodes_[v].label, v).second)
                throw std::invalid_argument("leaf labels must be nonempty and unique");
            begin_[v] = leaves_.size();
            leaves_.push_back(v);
            end_[v] = leaves_.size();
        }
        for (auto it = children.rbegin(); it != children.rend(); ++it) {
            const NodeId c = *it;
            if (c >= size() || seen[c])
                throw std::invalid_argument("edges must form a rooted tree");
            seen[c] = true;
            parents_[c] = v;
            depths_[c] = depths_[v] + 1;
            stack.push_back(c);
        }
    }
    if (preorder_.size() != size())
        throw std::invalid_argument("tree has unreachable nodes");
    for (auto it = preorder_.rbegin(); it != preorder_.rend(); ++it) {
        const auto& children = nodes_[*it].children;
        if (!children.empty()) {
            begin_[*it] = begin_[children.front()];
            end_[*it] = end_[children.back()];
        }
    }
}

NodeId Tree::leaf(const std::string& label) const {
    const auto found = labels_.find(label);
    if (found == labels_.end()) throw std::invalid_argument("unknown leaf label: " + label);
    return found->second;
}

void validate_leaf_sets(const Tree& first, const Tree& second) {
    if (first.leaf_count() != second.leaf_count())
        throw std::invalid_argument("trees have different leaf-label sets");
    for (NodeId v : first.leaves()) (void)second.leaf(first.node(v).label);
}

namespace {

// Iterative parser: even a long caterpillar uses no recursive C++ call stack.
class NewickParser {
public:
    explicit NewickParser(const std::string& text) : text_(text) {}

    Tree parse() {
        std::vector<Node> nodes;
        std::vector<NodeId> open;
        NodeId root = no_node;
        bool need_subtree = true;
        for (;;) {
            skip_space_and_comments();
            if (need_subtree) {
                if (peek() == '\0' || peek() == ')' || peek() == ',' || peek() == ';')
                    fail("expected a subtree");
                const NodeId v = nodes.size();
                nodes.push_back({{}, {}});
                if (open.empty()) {
                    if (root != no_node) fail("multiple roots");
                    root = v;
                } else {
                    nodes[open.back()].children.push_back(v);
                }
                if (peek() == '(') {
                    ++pos_;
                    open.push_back(v);
                    continue;
                }
                nodes[v].label = label();
                if (nodes[v].label.empty()) fail("missing leaf label");
                branch_length();
                need_subtree = false;
            }
            skip_space_and_comments();
            if (peek() == ',') {
                if (open.empty()) fail("comma outside an internal node");
                ++pos_;
                need_subtree = true;
            } else if (peek() == ')') {
                if (open.empty()) fail("unmatched closing parenthesis");
                ++pos_;
                const NodeId v = open.back();
                open.pop_back();
                nodes[v].label = label();
                branch_length();
            } else if (peek() == ';') {
                if (!open.empty()) fail("unclosed parenthesis");
                ++pos_;
                skip_space_and_comments();
                if (pos_ != text_.size()) fail("text after the tree");
                return Tree(std::move(nodes), root);
            } else {
                fail("expected comma, closing parenthesis or semicolon");
            }
        }
    }

private:
    const std::string& text_;
    std::size_t pos_ = 0;

    char peek() const { return pos_ == text_.size() ? '\0' : text_[pos_]; }
    [[noreturn]] void fail(const char* message) const {
        throw std::invalid_argument(std::string("Newick at byte ") +
                                    std::to_string(pos_) + ": " + message);
    }
    static bool space(char c) { return std::isspace(static_cast<unsigned char>(c)) != 0; }
    static bool delimiter(char c) {
        return c == '\0' || space(c) || c == '(' || c == ')' || c == ',' ||
               c == ':' || c == ';' || c == '[' || c == ']';
    }
    void skip_space_and_comments() {
        for (;;) {
            while (space(peek())) ++pos_;
            if (peek() != '[') return;
            std::size_t nesting = 0;
            do {
                const char c = peek();
                if (c == '\0') fail("unclosed comment");
                ++pos_;
                if (c == '[') ++nesting;
                if (c == ']') --nesting;
            } while (nesting);
        }
    }
    std::string label() {
        skip_space_and_comments();
        if (peek() != '\'') {
            const auto start = pos_;
            while (!delimiter(peek())) ++pos_;
            return text_.substr(start, pos_ - start);
        }
        ++pos_;
        std::string value;
        for (;;) {
            const char c = peek();
            if (c == '\0') fail("unclosed quoted label");
            ++pos_;
            if (c != '\'') {
                value.push_back(c);
            } else if (peek() == '\'') {
                value.push_back('\'');
                ++pos_;
            } else {
                return value;
            }
        }
    }
    void branch_length() {
        skip_space_and_comments();
        if (peek() != ':') return;
        ++pos_;
        skip_space_and_comments();
        const auto start = pos_;
        while (!delimiter(peek())) ++pos_;
        const auto value = text_.substr(start, pos_ - start);
        std::size_t consumed = 0;
        try {
            const double length = std::stod(value, &consumed);
            if (!std::isfinite(length) || consumed != value.size()) fail("invalid branch length");
        } catch (const std::exception&) {
            fail("invalid branch length");
        }
    }
};

}  // namespace

Tree Tree::from_newick(const std::string& text) { return NewickParser(text).parse(); }

}  // namespace soda2013
