/**
 * Master experiment runner — runs all experiments in sequence.
 * Compile:  cmake -B build && cmake --build build --config Release
 * Run:      ./build/run_all
 */
#include <cstdlib>
#include <iostream>

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  HBD Root-Finder — Complete Experimental Validation (C++/GMP)║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n\n";

    // Run each experiment as a separate process for clean measurement
    const char* exps[] = {
        "./build/exp8_baseline",
        "./build/exp4_speed",
        "./build/exp5_wilkinson",
        "./build/exp6_multiple",
        "./build/exp7_memory",
        "./build/exp9_vieta",
    };
    const char* names[] = {
        "Experiment 8: HBD vs HB-ED Baseline",
        "Experiment 4: Speed Quantification",
        "Experiment 5: Wilkinson W20",
        "Experiment 6: Multiple Roots",
        "Experiment 7: Memory Scaling",
        "Experiment 9: Pure Vieta vs Vieta+Cauchy",
    };

    int n_exps = sizeof(exps) / sizeof(exps[0]);
    for (int i = 0; i < n_exps; ++i) {
        std::cout << "\n  >> Running " << names[i] << "...\n";
        std::cout << "  -----------------------------------------------------------\n";
        int ret = std::system(exps[i]);
        if (ret != 0) {
            std::cerr << "  WARNING: " << names[i] << " returned " << ret << "\n";
        }
    }

    std::cout << "\n  All experiments complete.\n";
    return 0;
}
