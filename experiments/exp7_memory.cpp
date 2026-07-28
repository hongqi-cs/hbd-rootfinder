/**
 * Experiment 7: Memory Scaling
 *
 * Compares HBD O(n) vs Companion-QR O(n²) memory requirements.
 * n = 100, 500, 1000, 2000.
 */
#include "polynomial.h"
#include "hbd.h"
#include "timer.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>

#ifdef __linux__
#include <sys/resource.h>
#endif

using namespace hbd;

/// Get current RSS memory in KB
size_t peak_rss_kb() {
#ifdef __linux__
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            // Format: "VmRSS:    12345 kB"
            size_t pos = line.find(':');
            std::string val = line.substr(pos + 1);
            // Strip leading spaces and trailing " kB"
            val.erase(0, val.find_first_not_of(" \t"));
            val.erase(val.find(" kB"));
            return std::stoul(val);
        }
    }
#endif
    return 0;
}

int main() {
    std::cout << "=============================================================\n";
    std::cout << "  Experiment 7: Memory Scaling\n";
    std::cout << "  HBD O(n) vs Companion-QR O(n^2)\n";
    std::cout << "=============================================================\n\n";

    std::vector<int> ns = {100, 500, 1000, 2000};

    std::cout << std::left
              << std::setw(8)  << "n"
              << std::setw(14) << "HBD (KB)"
              << std::setw(14) << "HBD O(n)"
              << std::setw(14) << "QR O(n^2)"
              << "QR/HBD\n";
    std::cout << std::string(65, '-') << "\n";

    for (int n : ns) {
        set_gmp_precision(64);

        auto roots_d = random_roots(n, 5.0, 20240101 + n * 100);
        std::vector<mpf_class> roots_mp(n);
        for (int i = 0; i < n; ++i) roots_mp[i] = mpf_class(roots_d[i]);

        auto coeffs = poly_from_roots<mpf_class>(roots_mp);
        mpf_class bound(7.5);   // spread*1.5
        mpf_class tol("1e-8");

        size_t mem_before = peak_rss_kb();
        hbd_deflate_loop<mpf_class>(coeffs, bound, tol, std::min(400, n * 5), 100);
        size_t mem_after = peak_rss_kb();
        size_t delta = (mem_after > mem_before) ? (mem_after - mem_before) : 0;

        double hbd_est = double(n) * 8.0 / 1024.0;   // O(n) * 8 bytes -> KB
        double qr_est  = double(n) * double(n) * 8.0 / 1024.0;  // O(n^2) -> KB
        double ratio = (hbd_est > 0) ? qr_est / hbd_est : 0;

        std::cout << std::left
                  << std::setw(8)  << n
                  << std::setw(14) << std::fixed << std::setprecision(1) << double(delta)
                  << std::setw(14) << std::fixed << std::setprecision(1) << hbd_est
                  << std::setw(14) << std::fixed << std::setprecision(1) << qr_est
                  << std::fixed << std::setprecision(0) << ratio << "x\n";
    }

    std::cout << "\n  Paper Table 8.8:\n";
    std::cout << "    n=100:   4 KB   (HBD)  vs  80 KB   (QR)     = 20x\n";
    std::cout << "    n=500:   12 KB  (HBD)  vs  2 MB    (QR)     = 167x\n";
    std::cout << "    n=1000:  16 KB  (HBD)  vs  8 MB    (QR)     = 500x\n";
    std::cout << "    n=2000:  28 KB  (HBD)  vs  32 MB   (QR)     = 1143x\n";

    return 0;
}
