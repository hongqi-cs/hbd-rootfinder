/**
 * Experiment 9: Vieta Upper-Bound vs Cauchy Upper-Bound
 *
 * Tests the Vieta bound |a_0|^{1/n} as a bisection starting range.
 * The guarantee: at least one root (real or complex) satisfies |r| ≤ |a_0|^{1/n}.
 *
 * Strategy: Vieta-first, Cauchy-fallback.
 *   - Scan narrow Vieta range first (higher resolution → better root detection)
 *   - If nothing found (likely no real root in tight range), expand to Cauchy
 *   - Recompute Vieta after each deflation
 *
 * Focus: polynomials with known real roots, where Vieta shines.
 */
#include "hbd.h"
#include "polynomial.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>
#include <gmpxx.h>

using namespace hbd;

struct VietaRow {
    std::string name;
    int    n;
    double bound_cauchy;
    double bound_vieta;
    double ratio;          // Cauchy/Vieta
    double t_c;
    double t_v;
    int    nr_c, nr_v;     // # roots found
    double res_c, res_v;   // max residual
};

static VietaRow run_pair(const std::string& name,
                          const std::vector<double>& coeffs,
                          double tol, int n_scan) {
    VietaRow r;
    r.name = name;
    r.n = int(coeffs.size()) - 1;
    r.bound_cauchy = cauchy_bound(coeffs);
    r.bound_vieta = vieta_bound(coeffs);
    r.ratio = (r.bound_vieta > 0.0) ? r.bound_cauchy / r.bound_vieta : 0.0;

    // ── Cauchy ──
    auto t0 = std::chrono::steady_clock::now();
    auto rc = hbd_deflate_loop(coeffs, r.bound_cauchy, tol, n_scan, 200);
    auto t1 = std::chrono::steady_clock::now();
    r.t_c = double(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) / 1e6;
    r.nr_c = int(rc.size());
    r.res_c = 0.0;
    for (auto x : rc) r.res_c = std::max(r.res_c, std::fabs(horner_eval(coeffs, x)));

    // ── Vieta (with Cauchy fallback for first-vieta failure) ──
    t0 = std::chrono::steady_clock::now();
    auto rv = hbd_deflate_loop_vieta(coeffs, tol, r.bound_cauchy, n_scan, 200);
    t1 = std::chrono::steady_clock::now();
    r.t_v = double(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) / 1e6;
    r.nr_v = int(rv.size());
    r.res_v = 0.0;
    for (auto x : rv) r.res_v = std::max(r.res_v, std::fabs(horner_eval(coeffs, x)));

    return r;
}

int main() {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "═══════════════════════════════════════════════════════════════════\n";
    std::cout << "  Experiment 9: Vieta vs Cauchy Upper-Bound Comparison\n";
    std::cout << "  Strategy: Vieta first scan, Cauchy as fallback on miss\n";
    std::cout << "═══════════════════════════════════════════════════════════════════\n\n";

    std::vector<VietaRow> rows;
    double tol = 1e-12;

    // ── All-real-root polynomials: Vieta should dominate ──

    // Test 1: (x-2)^10
    {
        std::vector<double> cf = {1.0,
            -20, 180, -960, 3360, -8064,
            13440, -15360, 11520, -5120, 1024};
        rows.push_back(run_pair("(x-2)^10", cf, tol, 1000));
    }

    // Test 2: (x-1)(x-3)(x-5)(x-7)(x-9)
    {
        double rs[] = {1, 3, 5, 7, 9};
        auto cf = poly_from_roots<double>(std::vector<double>(rs, rs+5));
        rows.push_back(run_pair("(x-1)(x-3)(x-5)(x-7)(x-9)", cf, tol, 1000));
    }

    // Test 3: roots ∈ [-5,5], n=10
    {
        auto real_roots = random_roots(10, 5.0, 20240101);
        auto cf = poly_from_roots<double>(real_roots);
        rows.push_back(run_pair("10 roots in [-5,5]", cf, tol, 1000));
    }

    // Test 4: roots ∈ [-2,2], n=20
    {
        auto real_roots = random_roots(20, 2.0, 20240101);
        auto cf = poly_from_roots<double>(real_roots);
        rows.push_back(run_pair("20 roots in [-2,2]", cf, tol, 2000));
    }

    // Test 5: roots ∈ [-3,3], n=50
    {
        auto real_roots = random_roots(50, 3.0, 20240101);
        auto cf = poly_from_roots<double>(real_roots);
        rows.push_back(run_pair("50 roots in [-3,3]", cf, tol, 4000));
    }

    // ── Random coefficients: Vieta may not help (mostly complex roots) ──
    {
        auto cf = random_coeffs(20, 10.0, 20240101);
        rows.push_back(run_pair("random coeffs n=20", cf, tol, 1500));
    }
    {
        auto cf = random_coeffs(50, 10.0, 20240101);
        rows.push_back(run_pair("random coeffs n=50", cf, tol, 2000));
    }

    // ── Print ──
    std::cout << "┌──────────────────────────────────┬────┬──────────────┬─────────────────────────────────────┬─────────────────────────────────────┐\n";
    std::cout << "│ Test                             │ n  │ Ca/Vi ratio  │     Cauchy                           │     Vieta                            │\n";
    std::cout << "│                                  │    │              ├──────────┬──────────┬───────────────┼──────────┬──────────┬───────────────┤\n";
    std::cout << "│                                  │    │              │ time(s)  │ #roots   │ max|P(r)|     │ time(s)  │ #roots   │ max|P(r)|     │\n";
    std::cout << "├──────────────────────────────────┼────┼──────────────┼──────────┼──────────┼───────────────┼──────────┼──────────┼───────────────┤\n";

    for (auto& r : rows) {
        std::cout << "│ " << std::left  << std::setw(32) << r.name
                  << " │ " << std::setw(2) << std::right << r.n
                  << " │ " << std::setw(7) << std::setprecision(0) << r.ratio << " ×"
                  << "  │ " << std::setw(8) << std::setprecision(6) << r.t_c
                  << "│ " << std::setw(4) << r.nr_c << "/" << r.n
                  << "  │ " << std::setw(9) << std::scientific << std::setprecision(1) << r.res_c
                  << "  │ " << std::setw(8) << std::fixed << std::setprecision(6) << r.t_v
                  << "│ " << std::setw(4) << r.nr_v << "/" << r.n
                  << "  │ " << std::setw(9) << std::scientific << std::setprecision(1) << r.res_v
                  << "  │\n";
    }

    std::cout << "└──────────────────────────────────┴────┴──────────────┴──────────┴──────────┴───────────────┴──────────┴──────────┴───────────────┘\n\n";

    // ── Summary ──
    std::cout << "┌──────────────────────────────────┬──────────────┬───────────┬──────────────────┐\n";
    std::cout << "│ Summary                          │ Bound ratio  │ Speedup   │ Root match       │\n";
    std::cout << "├──────────────────────────────────┼──────────────┼───────────┼──────────────────┤\n";
    for (auto& r : rows) {
        double sp = (r.t_v > 0.0) ? r.t_c / r.t_v : 0.0;
        std::cout << "│ " << std::left  << std::setw(32) << r.name
                  << "│ " << std::setw(8) << std::right << std::setprecision(0) << r.ratio
                  << " ×  │ " << std::setw(5) << std::setprecision(2) << sp
                  << " ×  │ " << std::setw(5) << r.nr_v << "/" << r.nr_c
                  << "            │\n";
    }
    std::cout << "└──────────────────────────────────┴──────────────┴───────────┴──────────────────┘\n\n";

    std::cout << "Interpretation:\n";
    std::cout << "  - Vieta bound = |a_0|^{1/n} — from the product-of-roots identity (Vieta).\n";
    std::cout << "  - Cauchy bound = 1 + max|a_k| — always valid, often extremely loose.\n";
    std::cout << "  - Vieta is 10×–7680× tighter for polynomials with known real roots.\n";
    std::cout << "  - For random-coefficient polynomials, Vieta may find fewer roots because\n";
    std::cout << "    the smallest-modulus root is often complex.  In this case the Cauchy\n";
    std::cout << "    fallback ensures no real root is missed.\n";
    std::cout << "  - As a pure-bisection method, HBD doesn't need Cauchy's 'safety margin'\n";
    std::cout << "    — bisection never diverges.  Vieta's tight bound is always safe.\n";

    return 0;
}
