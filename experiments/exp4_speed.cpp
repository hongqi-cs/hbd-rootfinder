/**
 * Experiment 4: Speed Quantification
 *
 * HBD vs Newton (secant) vs Companion-QR at n=50, double precision.
 */
#include "polynomial.h"
#include "hbd.h"
#include "hbed.h"
#include "timer.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>

using namespace hbd;

// Simple Newton+secant: scan then refine with derivative estimation
std::vector<double> newton_secant(const std::vector<double>& coeffs,
                                   double bound, double tol) {
    std::vector<double> roots;
    auto curr = coeffs;
    while (curr.size() > 1) {
        int n_scan = 400;
        double L = -bound, R = bound;
        double step = (R - L) / n_scan;
        bool found = false;

        for (int i = 0; i < n_scan && !found; ++i) {
            double sL = L + step * i, sR = L + step * (i + 1);
            if (peval(curr, sL) * peval(curr, sR) > 0) continue;

            // Secant refinement from midpoint
            double x = (sL + sR) / 2.0;
            for (int it = 0; it < 100; ++it) {
                double fx = peval(curr, x);
                if (std::abs(fx) < tol) break;
                double df = (peval(curr, x + tol) - fx) / tol;
                if (std::abs(df) < 1e-40) break;
                double dx = fx / df;
                x -= dx;
                if (std::abs(dx) < tol) break;
            }

            roots.push_back(x);
            curr = synthetic_div(curr, x);
            found = true;
        }
        if (!found) break;
    }
    return roots;
}

int main() {
    std::cout << "=============================================================\n";
    std::cout << "  Experiment 4: Speed Quantification\n";
    std::cout << "  n=30, random int coeffs [-10,10], double precision\n";
    std::cout << "=============================================================\n\n";

    const int n = 30, n_runs = 10;
    std::vector<double> ht, nt;

    for (int run = 0; run < n_runs; ++run) {
        auto coeffs = random_coeffs(n, 10.0, 20240101 + run);
        double bound = 1.0;
        for (int i = 1; i <= n; ++i)
            bound = std::max(bound, 1.0 + std::abs(coeffs[i]));
        double tol = 1e-10;

        // HBD
        Timer t;
        hbd_deflate_loop(coeffs, bound, tol, 600, 200);
        ht.push_back(t.secs());

        // Newton+secant
        t.reset();
        newton_secant(coeffs, bound, tol);
        nt.push_back(t.secs());
    }

    double avg_h = std::accumulate(ht.begin(), ht.end(), 0.0) / n_runs;
    double avg_n = std::accumulate(nt.begin(), nt.end(), 0.0) / n_runs;

    std::cout << "  Avg time over " << n_runs << " runs:\n";
    std::cout << "    HBD (bisection):        " << std::fixed << std::setprecision(4)
              << avg_h << " s\n";
    std::cout << "    Newton+secant:          " << std::fixed << std::setprecision(4)
              << avg_n << " s\n";
    std::cout << "    Ratio (HBD/Newton):     " << std::fixed << std::setprecision(1)
              << (avg_n > 0 ? avg_h / avg_n : 0.0) << "x\n";

    std::cout << "\n  Paper Table 8.5 (n=50, double):\n";
    std::cout << "    HBD ~8.7x slower than Newton in standard double precision.\n";
    std::cout << "    (Note: our n=30 with GMP 64-bit may differ from paper's C++ double.)\n";

    return 0;
}
