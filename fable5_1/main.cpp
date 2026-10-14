// main.cpp — command line front end.
//
//   soda13 [--triplet|--quartet|--both] [-p P] [--naive] T1.nwk T2.nwk
//
// Prints the five cells of Fig. 3 and the (parameterized) distance.
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "distance.hpp"
#include "naive_distance.hpp"

using namespace soda13;

static std::string read_file(const char* path) {
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr, "cannot open %s\n", path);
        std::exit(2);
    }
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

static void print(const char* name, const CaseCounts& c, long double p) {
    // B + p (C + D) is a long double; print every significant digit the type
    // carries (LDBL_DIG, the decimal->long double->decimal round-trip width)
    // instead of a fixed number of decimals, so that a small p such as 0.01
    // on a single disagreeing quartet is not rounded away (trial 8,
    // Finding 2).  Integer-valued results still print as integers (%Lg
    // drops trailing zeros) up to 10^18; beyond that %Lg switches to
    // exponent notation, which is also where long double stops holding the
    // 128-bit totals exactly.
    std::printf("%s: A=%s B=%s C=%s D=%s E=%s distance=%s parameterized(p=%.*Lg)=%.*Lg\n", name,
                count_to_string(c.A).c_str(), count_to_string(c.B).c_str(), count_to_string(c.C).c_str(),
                count_to_string(c.D).c_str(), count_to_string(c.E).c_str(), count_to_string(c.distance()).c_str(),
                LDBL_DIG, p, LDBL_DIG, c.parameterized(p));
}

int main(int argc, char** argv) {
    bool trip = true, quart = true, naive = false;
    long double p = 1.0L;
    const char* files[2] = {nullptr, nullptr};
    int nf = 0;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--triplet")) {
            trip = true;
            quart = false;
        } else if (!std::strcmp(argv[i], "--quartet")) {
            trip = false;
            quart = true;
        } else if (!std::strcmp(argv[i], "--both")) {
            trip = quart = true;
        } else if (!std::strcmp(argv[i], "--naive")) {
            naive = true;
        } else if (!std::strcmp(argv[i], "-p") && i + 1 < argc) {
            p = std::strtold(argv[++i], nullptr);
        } else if (nf < 2) {
            files[nf++] = argv[i];
        } else {
            std::fprintf(stderr, "usage: %s [--triplet|--quartet|--both] [-p P] [--naive] T1.nwk T2.nwk\n", argv[0]);
            return 2;
        }
    }
    if (nf != 2) {
        std::fprintf(stderr, "usage: %s [--triplet|--quartet|--both] [-p P] [--naive] T1.nwk T2.nwk\n", argv[0]);
        return 2;
    }
    try {
        Tree t1 = Tree::from_newick(read_file(files[0]));
        Tree t2 = Tree::from_newick(read_file(files[1]));
        if (naive) assign_leaf_ids(t1, t2);  // the oracles work on leaf ids directly
        if (trip) print(naive ? "triplets(naive)" : "triplets", naive ? naive_triplet_counts(t1, t2) : triplet_case_counts(t1, t2), p);
        if (quart) print(naive ? "quartets(naive)" : "quartets", naive ? naive_quartet_counts(t1, t2) : quartet_case_counts(t1, t2), p);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
