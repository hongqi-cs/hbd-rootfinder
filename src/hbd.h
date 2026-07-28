#pragma once
#include "polynomial.h"
#include <vector>
#include <utility>
#include <algorithm>

namespace hbd {

// ════════════════════════════════════════════════════════════════
//  HBD (Horner-Bisection Deflation)
//
//  Uses Phi(x) for sign-change detection and bisection.
//  When bisection converges, the Horner intermediate sequence
//  [f_1, ..., f_{n-1}] directly gives the deflated quotient.
//  Root r = -x* (where x* is the Phi root).
//
//  NOTE: All T(int) replaced with T(double) for GMP compatibility.
// ════════════════════════════════════════════════════════════════

/// Bisection on Phi(x) within [L, R]
/// Returns (x*, quotient_f_seq, success)
template <typename T>
std::tuple<T, std::vector<T>, bool>
bisect_hbd(const std::vector<T>& coeffs, T L, T R,
           T tol, int max_iter) {
    auto [pL, fL] = phi_eval(coeffs, L);
    auto [pR, fR] = phi_eval(coeffs, R);

    if (pL * pR > T(0.0))
        return {T(0.0), {}, false};

    if (myabs(pL) < tol)
        return {L, std::move(fL), true};
    if (myabs(pR) < tol)
        return {R, std::move(fR), true};

    for (int iter = 0; iter < max_iter; ++iter) {
        T mid = (L + R) / T(2.0);
        auto [pm, fm] = phi_eval(coeffs, mid);

        if (myabs(pm) < tol)
            return {mid, std::move(fm), true};

        if (pL * pm < T(0.0)) {
            R = mid; pR = pm;
        } else {
            L = mid; pL = pm;
        }

        if (myabs(T(R - L)) < tol) {
            T mid2 = (L + R) / T(2.0);
            auto [_, fm2] = phi_eval(coeffs, mid2);
            return {mid2, std::move(fm2), true};
        }
    }

    T mid = (L + R) / T(2.0);
    auto [phi_mid, fm] = phi_eval(coeffs, mid);
    bool ok = myabs(phi_mid) < tol * T(10.0);
    return {mid, std::move(fm), ok};
}

/// Scan + bisect + deflate loop using HBD
/// coeffs: monic [a_n=1, ..., a_0]
/// bound: root search range [-bound, bound]
template <typename T>
std::vector<T> hbd_deflate_loop(const std::vector<T>& coeffs_init,
                                 T bound, T tol,
                                 int n_scan = 200, int max_iter = 200) {
    std::vector<T> roots;
    std::vector<T> curr = coeffs_init;  // copy

    while (curr.size() > 1) {
        T L = -bound, R = bound;
        T step = (R - L) / T(double(n_scan));
        bool found = false;

        // Phase 1: scan for sign changes in Phi(x)
        T prev_x = L;
        auto [prev_phi, _] = phi_eval(curr, prev_x);
        T prev_sign = (myabs(prev_phi) < tol) ? T(0.0) :
                ((prev_phi > T(0.0)) ? T(1.0) : T(-1.0));

        for (int i = 1; i <= n_scan; ++i) {
            T curr_x = L + step * T(double(i));
            auto [curr_phi, _] = phi_eval(curr, curr_x);
            T curr_sign = (myabs(curr_phi) < tol) ? T(0.0) :
                    ((curr_phi > T(0.0)) ? T(1.0) : T(-1.0));

            if (prev_sign * curr_sign < T(0.0) ||
                (curr_sign == T(0.0) && prev_sign != T(0.0))) {
                auto [x_star, fseq, ok] = bisect_hbd(
                    curr, prev_x, curr_x, tol, max_iter);
                if (ok && !fseq.empty()) {
                    roots.push_back(-x_star);     // r = -x*
                    curr = std::move(fseq);        // deflate
                    found = true;
                    break;
                }
            }
            prev_x = curr_x;
            prev_sign = curr_sign;
        }

        // Phase 2: finer fallback scan
        if (!found) {
            int n_fine = n_scan * 4;
            T step2 = (R - L) / T(double(n_fine));
            for (int i = 1; i <= n_fine; ++i) {
                T sL = L + step2 * T(double(i - 1));
                T sR = L + step2 * T(double(i));
                auto [pL, _] = phi_eval(curr, sL);
                auto [pR, __] = phi_eval(curr, sR);
                if (pL * pR > T(0.0)) continue;

                auto [x_star, fseq, ok] = bisect_hbd(
                    curr, sL, sR, tol, max_iter);
                if (ok && !fseq.empty()) {
                    roots.push_back(-x_star);
                    curr = std::move(fseq);
                    found = true;
                    break;
                }
            }
        }

        if (!found) break;  // no more real roots
    }

    return roots;
}

/// Vieta-optimised variant: bound is recomputed from |a_0|^{1/n}
/// after each deflation, progressively tightening the search range.
///
/// Strategy (Vieta-first, Cauchy-fallback):
///   1. Compute Vieta bound R_V = |a_0|^{1/k} (k = current degree)
///   2. Scan [-R_V, R_V] at high resolution
///   3. If nothing found, expand to Cauchy fallback (full range)
///   4. After each deflation, recompute R_V for the reduced polynomial
///
/// This ensures tight range for well-behaved real-root polynomials
/// while guaranteeing all real roots are found via Cauchy fallback.
template <typename T>
std::vector<T> hbd_deflate_loop_vieta(const std::vector<T>& coeffs_init,
                                       T tol, T fallback_bound,
                                       int n_scan = 200, int max_iter = 200) {
    std::vector<T> roots;
    std::vector<T> curr = coeffs_init;

    while (curr.size() > 1) {
        int k = int(curr.size()) - 1;
        T a0_abs = myabs(curr.back());
        T bound;
        bool use_vieta = true;

        if (a0_abs == T(0.0)) {
            roots.push_back(T(0.0));
            curr.pop_back();
            continue;
        } else {
            bound = vieta_bound(curr);
            if (bound > fallback_bound) {
                bound = fallback_bound;
                use_vieta = false;
            }
        }

        // Helper lambda: scan [L, R] with given n_pts
        auto scan_and_bisect = [&](T L, T R, int n_pts, bool& found) {
            T step = (R - L) / T(double(n_pts));
            T prev_x = L;
            auto [prev_phi, _] = phi_eval(curr, prev_x);
            T prev_sign = (myabs(prev_phi) < tol) ? T(0.0) :
                    ((prev_phi > T(0.0)) ? T(1.0) : T(-1.0));

            for (int i = 1; i <= n_pts; ++i) {
                T x = L + step * T(double(i));
                auto [phi, f] = phi_eval(curr, x);
                T sign = (myabs(phi) < tol) ? T(0.0) :
                        ((phi > T(0.0)) ? T(1.0) : T(-1.0));
                if (prev_sign * sign < T(0.0) ||
                    (sign == T(0.0) && prev_sign != T(0.0))) {
                    auto [x_star, fseq, ok] = bisect_hbd(
                        curr, prev_x, x, tol, max_iter);
                    if (ok && !fseq.empty()) {
                        roots.push_back(-x_star);
                        curr = std::move(fseq);
                        found = true;
                        return;
                    }
                }
                prev_x = x;
                prev_sign = sign;
            }
        };

        bool found = false;

        // ── Step 1: Vieta-range scan (tight, high resolution) ──
        scan_and_bisect(-bound, bound, n_scan, found);

        // ── Step 2: fine scan within Vieta range (×4 resolution) ──
        if (!found && use_vieta) {
            scan_and_bisect(-bound, bound, n_scan * 4, found);
        }

        // ── Step 3: Cauchy fallback — expand to full safe range ──
        if (!found && use_vieta && bound < fallback_bound) {
            scan_and_bisect(-fallback_bound, fallback_bound, n_scan, found);
            if (!found) {
                scan_and_bisect(-fallback_bound, fallback_bound, n_scan * 4, found);
            }
        }

        if (!found) break;
    }

    return roots;
}

} // namespace hbd
