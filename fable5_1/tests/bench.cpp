// tests/bench.cpp — running-time scaling check (not an oracle test).
//
//   bench [mode] [nmax] [shape]
//     mode  : t (triplets), q (quartets), tq
//     shape : b (binary-ish random), m (mixed degrees), s (star-ish, higher degree),
//             S (root star of sqrt(n) random binary subtrees: d ~ sqrt(n)),
//             c (caterpillar: unbalanced, stresses the HDT balance)
//
// Prints n, max degree, wall time and time / (n log n) so that the
// O(n log n) (triplets) and O(d n log n) (quartets) bounds of the paper
// can be eyeballed.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../distance.hpp"
#include "test_util.hpp"

using namespace soda13;

int main(int argc, char** argv) {
    const char* mode = argc > 1 ? argv[1] : "tq";
    int nmax = argc > 2 ? std::atoi(argv[2]) : 100000;
    char shape = argc > 3 ? argv[3][0] : 'b';
    std::mt19937 rng(99);
    for (int n = 1000; n <= nmax; n *= 2) {
        double bias = shape == 'b' ? 1.0 : shape == 'm' ? 0.5 : 0.05;
        Tree t1, t2;
        if (shape == 'S') {
            t1 = test::star_of_binaries(n, rng);
            t2 = test::random_tree(n, rng, 1.0);
        } else if (shape == 'c') {
            t1 = test::caterpillar(n);
            t2 = test::caterpillar(n);
        } else {
            t1 = test::random_tree(n, rng, bias);
            t2 = test::random_tree(n, rng, bias);
        }
        test::shuffle_labels(t2, rng);
        assign_leaf_ids(t1, t2);
        t1.contract_unary();
        t2.contract_unary();
        const int d = std::max(t1.degree_max(), t2.degree_max());
        const double nlogn = n * std::log2((double)n);
        if (std::strchr(mode, 't')) {
            auto a = std::chrono::steady_clock::now();
            CaseCounts c = triplet_case_counts(t1, t2);
            auto b = std::chrono::steady_clock::now();
            double s = std::chrono::duration<double>(b - a).count();
            std::printf("triplet n=%7d d=%4d time=%8.3fs  time/(n log n)=%.3e  dist=%s\n", n, d, s, s / nlogn,
                        count_to_string(c.distance()).c_str());
        }
        if (std::strchr(mode, 'q')) {
            auto a = std::chrono::steady_clock::now();
            CaseCounts c = quartet_case_counts(t1, t2);
            auto b = std::chrono::steady_clock::now();
            double s = std::chrono::duration<double>(b - a).count();
            std::printf("quartet n=%7d d=%4d time=%8.3fs  time/(n log n)=%.3e  time/(d n log n)=%.3e  dist=%s\n", n, d,
                        s, s / nlogn, s / (d * nlogn), count_to_string(c.distance()).c_str());
        }
        std::fflush(stdout);
    }
    return 0;
}
