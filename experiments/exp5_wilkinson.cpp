/**
 * Experiment 5: Wilkinson W20
 *
 * W20(x) = prod_{i=1}^{20} (x - i)
 * Tests HBD's robustness on extremely ill-conditioned polynomials.
 * Uses GMP 256-bit precision.
 */
#include "polynomial.h"
#include "hbd.h"
#include "timer.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace hbd;

int main() {
    std::cout << "=============================================================\n";
    std::cout << "  Experiment 5: Wilkinson W20\n";
    std::cout << "  W20(x) = prod_{i=1}^{20} (x-i), GMP 256-bit\n";
    std::cout << "=============================================================\n\n";

    set_gmp_precision(256);

    // Build W20 from roots
    std::vector<mpf_class> true_roots(20);
    for (int i = 0; i < 20; ++i) true_roots[i] = mpf_class(double(i + 1));
    auto coeffs = poly_from_roots<mpf_class>(true_roots);

    mpf_class bound(21.0);       // roots in [1, 20]
    mpf_class tol("1e-20");
    int n_scan = 2000;

    // HBD
    Timer t;
    auto roots_hbd = hbd_deflate_loop<mpf_class>(coeffs, bound, tol, n_scan, 400);
    double t_hbd = t.secs();

    // Match to true roots
    std::vector<double> dists;
    for (int tn = 1; tn <= 20; ++tn) {
        mpf_class target = mpf_class(double(tn));
        double best = 1e99;
        for (const auto& r : roots_hbd) {
            double d = std::abs(mpf_class(r - target).get_d());
            if (d < 0.5 && d < best) best = d;
        }
        if (best < 1e99) dists.push_back(best);
    }

    mpf_class res = max_residual<mpf_class>(coeffs, roots_hbd);
    double max_dist = dists.empty() ? 1e99 : *std::max_element(dists.begin(), dists.end());

    std::cout << "  HBD (256-bit):\n";
    std::cout << "    Roots found: " << roots_hbd.size() << "/20\n";
    std::cout << "    Max |root - true|: " << std::scientific << max_dist << "\n";
    std::cout << "    Max residual:      " << std::scientific << res.get_d() << "\n";
    std::cout << "    Time:              " << std::fixed << std::setprecision(4) << t_hbd << " s\n";

    std::cout << "\n  Paper Table 8.6:\n";
    std::cout << "    HBD (256-bit):        20/20, err=1.8e-5,  92.7 ms\n";
    std::cout << "    Newton (double):      0/20  (diverges)\n";
    std::cout << "    Companion-QR (double): 6/20,  err=1.4e-4\n";
    std::cout << "    MPSolve (256-bit):    20/20, err=1.2e-16, 0.3 ms\n";

    return 0;
}
