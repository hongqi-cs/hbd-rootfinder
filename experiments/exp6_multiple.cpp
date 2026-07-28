/**
 * Experiment 6: Multiple Roots
 *
 * P(x) = (x-1)^3 (x+2)  vs  P(x) = (x-1)^5 (x+2)^3
 * Tests HBD convergence on multiple roots.
 * GMP 100-bit.
 */
#include "polynomial.h"
#include "hbd.h"
#include "timer.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace hbd;

void test_poly(const std::vector<mpf_class>& roots_true, const std::string& desc) {
    set_gmp_precision(100);
    auto coeffs = poly_from_roots<mpf_class>(roots_true);

    mpf_class bound(4.0);
    mpf_class tol("1e-12");
    int n_scan = 1200;

    // HBD with convergence tracking
    auto curr = coeffs;
    std::vector<mpf_class> roots_hbd;
    std::vector<int> conv_steps;

    while (curr.size() > 1) {
        int n = int(curr.size()) - 1;
        mpf_class L = -bound, R = bound;
        mpf_class step = (R - L) / mpf_class(n_scan);
        bool found = false;

        for (int i = 0; i < n_scan; ++i) {
            mpf_class sL = L + step * mpf_class(double(i));
            mpf_class sR = L + step * mpf_class(double(i + 1));
            auto [pL, _] = phi_eval<mpf_class>(curr, sL);
            auto [pR, __] = phi_eval<mpf_class>(curr, sR);
            if (pL * pR > mpf_class(0.0)) continue;

            int bisect_iters = 0;
            mpf_class mid;
            for (int iter = 0; iter < 400; ++iter) {
                ++bisect_iters;
                mid = (sL + sR) / mpf_class(2.0);
                auto [pm, fm] = phi_eval<mpf_class>(curr, mid);
                if (abs(pm) < tol) {
                    roots_hbd.push_back(-mid);
                    curr = std::move(fm);
                    conv_steps.push_back(bisect_iters);
                    found = true; break;
                }
                if (pL * pm < mpf_class(0.0)) { sR = mid; pR = pm; }
                else { sL = mid; pL = pm; }
                if (abs(sR - sL) < tol) {
                    mid = (sL + sR) / mpf_class(2.0);
                    auto [_, fm] = phi_eval<mpf_class>(curr, mid);
                    roots_hbd.push_back(-mid);
                    curr = std::move(fm);
                    conv_steps.push_back(bisect_iters);
                    found = true; break;
                }
            }
            if (found) break;
        }
        if (!found) break;
    }

    mpf_class res = max_residual<mpf_class>(coeffs, roots_hbd);

    // Cluster near 1 and -2
    std::vector<double> cluster_1, cluster_n2;
    for (const auto& r : roots_hbd) {
        double rd = r.get_d();
        if (std::abs(rd - 1.0) < 0.5) cluster_1.push_back(rd);
        if (std::abs(rd + 2.0) < 0.5) cluster_n2.push_back(rd);
    }

    double avg_iters = 0;
    if (!conv_steps.empty())
        avg_iters = std::accumulate(conv_steps.begin(), conv_steps.end(), 0.0) / conv_steps.size();

    std::cout << "  " << desc << ":\n";
    std::cout << "    Roots found: " << roots_hbd.size()
              << " (" << cluster_1.size() << " near 1, "
              << cluster_n2.size() << " near -2)\n";
    std::cout << "    Residual:    " << std::scientific << res.get_d() << "\n";

    if (!cluster_1.empty()) {
        std::sort(cluster_1.begin(), cluster_1.end());
        double spread = cluster_1.back() - cluster_1.front();
        double avg_err = 0;
        for (double x : cluster_1) avg_err += std::abs(x - 1.0);
        avg_err /= cluster_1.size();
        std::cout << "    Cluster @1:  spread=" << std::scientific << spread
                  << ", avg error=" << avg_err << "\n";
    }
    if (!cluster_n2.empty()) {
        std::sort(cluster_n2.begin(), cluster_n2.end());
        double spread = cluster_n2.back() - cluster_n2.front();
        double avg_err = 0;
        for (double x : cluster_n2) avg_err += std::abs(x + 2.0);
        avg_err /= cluster_n2.size();
        std::cout << "    Cluster @-2: spread=" << std::scientific << spread
                  << ", avg error=" << avg_err << "\n";
    }
    std::cout << "    Avg bisection iters/root: " << std::fixed
              << std::setprecision(1) << avg_iters << " (linear, ~1 bit/step)\n\n";
}

int main() {
    std::cout << "=============================================================\n";
    std::cout << "  Experiment 6: Multiple Roots\n";
    std::cout << "  HBD on (x-1)^k (x+2)^m, GMP 100-bit\n";
    std::cout << "=============================================================\n\n";

    // Test 1: (x-1)^3 (x+2)
    std::vector<mpf_class> roots1 = {
        mpf_class(1.0), mpf_class(1.0), mpf_class(1.0), mpf_class(-2.0)
    };
    test_poly(roots1, "P(x) = (x-1)^3 (x+2)");

    // Test 2: (x-1)^5 (x+2)^3 (from paper)
    std::vector<mpf_class> roots2 = {
        mpf_class(1.0), mpf_class(1.0), mpf_class(1.0), mpf_class(1.0), mpf_class(1.0),
        mpf_class(-2.0), mpf_class(-2.0), mpf_class(-2.0)
    };
    test_poly(roots2, "P(x) = (x-1)^5 (x+2)^3");

    std::cout << "  Paper Table 8.7:\n";
    std::cout << "    Newton:   linear conv, ~1-1/5=0.8, error ~1e-5\n";
    std::cout << "    HBD:      linear conv (bisection), 1 bit/step, error ~1e-12\n";

    return 0;
}
