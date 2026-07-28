#pragma once
#include "polynomial.h"
#include <vector>
#include <tuple>

namespace hbd {

// ════════════════════════════════════════════════════════════════
//  HB-ED (Horner-Bisection + Explicit Deflation) — baseline
//
//  Standard approach:
//  1. Bisection on P(t) to locate root t0
//  2. After convergence, run an additional synthetic division
//     to get the quotient (extra O(n) Horner pass)
//
//  This is the baseline against which HBD's single-stage
//  error injection advantage is measured.
//
//  NOTE: All T(int) replaced with T(double) for GMP compatibility.
// ════════════════════════════════════════════════════════════════

/// Bisection on P(t) within [L, R]; does NOT return quotient
template <typename T>
std::tuple<T, bool>
bisect_hbed(const std::vector<T>& coeffs, T L, T R,
            T tol, int max_iter) {
    T pL = peval(coeffs, L);
    T pR = peval(coeffs, R);

    if (pL * pR > T(0.0))
        return {T(0.0), false};

    if (myabs(pL) < tol)
        return {L, true};
    if (myabs(pR) < tol)
        return {R, true};

    T pL_val = pL;
    for (int iter = 0; iter < max_iter; ++iter) {
        T mid = (L + R) / T(2.0);
        T pm = peval(coeffs, mid);

        if (myabs(pm) < tol)
            return {mid, true};

        if (pL_val * pm < T(0.0)) {
            R = mid; pR = pm;
        } else {
            L = mid; pL_val = pm;
        }

        if (myabs(T(R - L)) < tol) {
            return {(L + R) / T(2.0), true};
        }
    }

    return {(L + R) / T(2.0), true};
}

/// Scan + bisect + explicit deflate loop using HB-ED
template <typename T>
std::vector<T> hbed_deflate_loop(const std::vector<T>& coeffs_init,
                                  T bound, T tol,
                                  int n_scan = 200, int max_iter = 200) {
    std::vector<T> roots;
    std::vector<T> curr = coeffs_init;

    while (curr.size() > 1) {
        T L = -bound, R = bound;
        T step = (R - L) / T(double(n_scan));
        bool found = false;

        // Phase 1: scan for sign changes in P(t)
        T prev_x = L;
        T prev_val = peval(curr, prev_x);
        T prev_sign = (myabs(prev_val) < tol) ? T(0.0) :
                ((prev_val > T(0.0)) ? T(1.0) : T(-1.0));

        for (int i = 1; i <= n_scan; ++i) {
            T curr_x = L + step * T(double(i));
            T curr_val = peval(curr, curr_x);
            T curr_sign = (myabs(curr_val) < tol) ? T(0.0) :
                    ((curr_val > T(0.0)) ? T(1.0) : T(-1.0));

            if (prev_sign * curr_sign < T(0.0) ||
                (curr_sign == T(0.0) && prev_sign != T(0.0))) {
                // Bisection on P(t)
                T subL = prev_x, subR = curr_x;
                T psubL = prev_val;
                bool bfound = false;
                T t0;

                for (int iter = 0; iter < max_iter; ++iter) {
                    T mid = (subL + subR) / T(2.0);
                    T pm = peval(curr, mid);

                    if (myabs(pm) < tol) {
                        t0 = mid; bfound = true; break;
                    }
                    if (psubL * pm < T(0.0)) {
                        subR = mid;
                    } else {
                        subL = mid; psubL = pm;
                    }
                    if (myabs(T(subR - subL)) < tol) {
                        t0 = (subL + subR) / T(2.0);
                        bfound = true; break;
                    }
                }
                if (!bfound)
                    t0 = (subL + subR) / T(2.0);

                roots.push_back(t0);
                // Explicit deflation — EXTRA O(n) step
                curr = synthetic_div(curr, t0);
                found = true;
                break;
            }
            prev_x = curr_x;
            prev_sign = curr_sign;
            prev_val = curr_val;
        }

        // Phase 2: finer fallback scan
        if (!found) {
            int n_fine = n_scan * 4;
            T step2 = (R - L) / T(double(n_fine));
            for (int i = 1; i <= n_fine; ++i) {
                T sL = L + step2 * T(double(i - 1));
                T sR = L + step2 * T(double(i));
                if (peval(curr, sL) * peval(curr, sR) > T(0.0))
                    continue;
                auto [t0, ok] = bisect_hbed(curr, sL, sR, tol, max_iter);
                if (ok) {
                    roots.push_back(t0);
                    curr = synthetic_div(curr, t0);
                    found = true;
                    break;
                }
            }
        }

        if (!found) break;
    }

    return roots;
}

} // namespace hbd
