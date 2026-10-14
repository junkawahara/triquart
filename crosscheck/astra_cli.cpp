// astra_cli.cpp -- minimal command-line front end for astra6, used only by
// crosscheck/. Prints the five counts A..E in the same layout as
// fable5_1/build/soda13 (without the distance columns).
//
//   astra_cli T1.nwk T2.nwk
#include <fstream>
#include <iostream>
#include <sstream>

#include "quartet_distance.hpp"
#include "triplet_distance.hpp"

using namespace soda2013;

static std::string read_file(const char* path) {
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

static void print(const char* name, const DistanceCounts& c) {
    std::cout << name << ": A=" << c.A << " B=" << c.B << " C=" << c.C << " D=" << c.D << " E=" << c.E << "\n";
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " T1.nwk T2.nwk\n";
        return 2;
    }
    try {
        const auto t1 = Tree::from_newick(read_file(argv[1]));
        const auto t2 = Tree::from_newick(read_file(argv[2]));
        print("triplets", triplet_counts(t1, t2));
        print("quartets", quartet_counts(t1, t2));
    } catch (const std::exception& e) {
        std::cout << "ERROR " << e.what() << "\n";
        return 1;
    }
}
