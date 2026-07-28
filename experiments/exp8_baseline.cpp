/**
 * Experiment 8: HBD vs HB-ED — Baseline Comparison (CORE)
 *
 * Reproduces Table 8.9 from the paper.
 * Two test modes:
 *   Small n (≤20): poly_from_roots (all roots real, verify correctness)
 *   Large n (≥50): random integer coeffs (Cauchy bound stays ~11, test
 *                  speed + residual on real roots found)
 */

#include "polynomial.h"
#include "hbd.h"
#include "hbed.h"
#include "timer.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace hbd;

struct Exp8Row {
    int n;
    int nh_hbd, nh_hbed;
    double t_hbd, t_hbed;
    double res_hbd, res_hbed;
    double speedup, prec_ratio;
};

// ─── Small-n test: poly_from_roots, all roots known real ─────────
Exp8Row run_exp8_small(int n, int n_scan, double tol_d, uint64_t seed) {
    auto roots_true = random_roots(n, 2.0, seed);
    auto coeffs = poly_from_roots<double>(roots_true);
    double bound = std::max(cauchy_bound<double>(coeffs), 0.0);
    if (bound > 10.0) bound = 4.0;  // roots ∈ [-2,2]

    Timer t;
    auto roots_hbd = hbd_deflate_loop(coeffs, bound, tol_d, n_scan, 200);
    double th = t.secs();

    t.reset();
    auto roots_hbed = hbed_deflate_loop(coeffs, bound, tol_d, n_scan, 200);
    double te = t.secs();

    double rh = max_residual(coeffs, roots_hbd);
    double re = max_residual(coeffs, roots_hbed);

    int nh = int(roots_hbd.size()), ne = int(roots_hbed.size());
    double sp = (th > 0) ? te / th : 0.0;
    double pr = (rh > 1e-300 && re > 1e-300) ? re / rh : 0.0;

    return {n, nh, ne, th, te, rh, re, sp, pr};
}

// ─── Large-n test: random integer coeffs, Cauchy bound stays small ───
Exp8Row run_exp8_large(int n, int n_scan, double tol_d, uint64_t seed) {
    auto coeffs = random_coeffs(n, 10.0, seed);
    double bound = cauchy_bound<double>(coeffs);

    Timer t;
    auto roots_hbd = hbd_deflate_loop(coeffs, bound, tol_d, n_scan, 200);
    double th = t.secs();

    t.reset();
    auto roots_hbed = hbed_deflate_loop(coeffs, bound, tol_d, n_scan, 200);
    double te = t.secs();

    double rh = max_residual(coeffs, roots_hbd);
    double re = max_residual(coeffs, roots_hbed);

    int nh = int(roots_hbd.size()), ne = int(roots_hbed.size());
    double sp = (th > 0) ? te / th : 0.0;
    double pr = (rh > 1e-300 && re > 1e-300) ? re / rh : 0.0;

    return {n, nh, ne, th, te, rh, re, sp, pr};
}

// ─── Main ──────────────────────────────────────────────────────────
int main() {
    std::cout << "============================================================================\n";
    std::cout << "  Experiment 8: HBD vs HB-ED Baseline Comparison\n";
    std::cout << "  Small n (≤20): poly_from_roots, verify correctness\n";
    std::cout << "  Large n (≥50): random integer coeffs [-10,10], speed + residual\n";
    std::cout << "============================================================================\n\n";

    struct Cfg { int n, n_scan; double tol; uint64_t seed; bool small; };
    std::vector<Cfg> configs = {
        {10,   200,  1e-12, 20240101, true},
        {20,   400,  1e-10, 20240102, true},
        {50,   400,  1e-10, 20240103, false},
        {100,  800,  1e-10, 20240104, false},
        {200,  1200, 1e-9,  20240105, false},
        {500,  1600, 1e-8,  20240106, false},
    };

    std::vector<Exp8Row> rows;

    for (const auto& cfg : configs) {
        std::cout << "  n=" << cfg.n << " (scan=" << cfg.n_scan
                  << ", tol=" << std::scientific << cfg.tol << ")..."
                  << std::flush;

        Exp8Row row;
        if (cfg.small) {
            row = run_exp8_small(cfg.n, cfg.n_scan, cfg.tol, cfg.seed);
        } else {
            row = run_exp8_large(cfg.n, cfg.n_scan, cfg.tol, cfg.seed);
        }

        std::cout << " done\n";
        std::cout << "    HBD:   " << row.nh_hbd << "/" << row.n
                  << " roots, " << row.t_hbd << "s, residual="
                  << std::scientific << row.res_hbd << "\n";
        std::cout << "    HB-ED: " << row.nh_hbed << "/" << row.n
                  << " roots, " << row.t_hbed << "s, residual="
                  << row.res_hbed << "\n";
        if (row.prec_ratio > 0)
            std::cout << "    Speed: " << std::fixed << std::setprecision(2)
                      << row.speedup << "x, Precision: HBD " << std::fixed
                      << std::setprecision(1) << row.prec_ratio
                      << "x more precise\n";
        else
            std::cout << "    Speed: " << std::fixed << std::setprecision(2)
                      << row.speedup << "x, Precision: N/A\n";

        rows.push_back(row);
    }

    // ── Summary table ──
    std::cout << "\n" << std::string(100, '=') << "\n";
    std::cout << "  Summary Table (cf. Paper Table 8.9)\n";
    std::cout << std::string(100, '=') << "\n";
    std::cout << std::left
              << std::setw(6)  << "n"
              << std::setw(9)  << "HBD(s)"
              << std::setw(9)  << "HB-ED(s)"
              << std::setw(12) << "HBD err"
              << std::setw(12) << "HB-ED err"
              << std::setw(8)  << "Speed"
              << std::setw(12) << "Prec-ratio"
              << std::setw(14) << "HBD found"
              << "HB-ED found\n";
    std::cout << std::string(100, '-') << "\n";

    for (const auto& r : rows) {
        std::cout << std::left
                  << std::setw(6)  << r.n
                  << std::setw(9)  << std::fixed << std::setprecision(3) << r.t_hbd
                  << std::setw(9)  << std::fixed << std::setprecision(3) << r.t_hbed
                  << std::setw(12) << std::scientific << std::setprecision(2) << r.res_hbd
                  << std::setw(12) << std::scientific << std::setprecision(2) << r.res_hbed
                  << std::setw(8)  << std::fixed << std::setprecision(2) << r.speedup << "x"
                  << std::setw(12);
        if (r.prec_ratio > 0)
            std::cout << std::fixed << std::setprecision(1) << r.prec_ratio << "x";
        else
            std::cout << "N/A";
        std::cout << std::setw(14) << (std::to_string(r.nh_hbd) + "/" + std::to_string(r.n))
                  << (std::to_string(r.nh_hbed) + "/" + std::to_string(r.n)) << "\n";
    }

    std::cout << "\n  Paper Table 8.9 (C++, i9-13900H, GCC 13.2 -O3):\n";
    std::cout << "    n=10:   0.002s / 0.003s,   4.2e-16 / 5.1e-16,  1.5x speed\n";
    std::cout << "    n=50:   0.048s / 0.072s,   8.7e-16 / 2.3e-15,  1.5x speed\n";
    std::cout << "    n=100:  0.195s / 0.291s,   3.1e-15 / 9.8e-14,  1.5x speed, 31x precise\n";
    std::cout << "    n=500:  4.87s  / 7.36s,    2.1e-14 / 1.4e-11,  1.5x speed, 670x precise\n";

    return 0;
}
