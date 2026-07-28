#pragma once
#include <vector>
#include <random>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <gmpxx.h>

namespace hbd {

// myabs: two overloads — no template to avoid GMP expression template deduction
inline double myabs(double x) { return std::fabs(x); }
inline mpf_class myabs(mpf_class x) { return abs(x); }

// ════════════════════════════════════════════════════════════════
//  Polynomial generation
// ════════════════════════════════════════════════════════════════

/// Build monic polynomial coeffs [a_n=1, ..., a_0] from roots
template <typename T>
std::vector<T> poly_from_roots(const std::vector<T>& roots) {
    size_t n = roots.size();
    std::vector<T> coeffs(n + 1, T(0.0));
    coeffs[0] = T(1.0);
    for (const auto& r : roots) {
        std::vector<T> new_c(n + 1, T(0.0));
        for (size_t i = 0; i <= n; ++i) {
            T prev = (i > 0) ? coeffs[i - 1] : T(0.0);
            new_c[i] = coeffs[i] - r * prev;
        }
        coeffs = std::move(new_c);
    }
    return coeffs;
}

/// Random integer coefficients in [-M, M] (non-zero), with fixed seed
inline std::vector<double> random_coeffs(int n, double M = 10.0, uint64_t seed = 20240101) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> dist(-int(M), int(M));
    std::vector<double> coeffs(n + 1);
    coeffs[0] = 1.0;
    for (int i = 1; i <= n; ++i) {
        int c = dist(rng);
        while (c == 0) c = dist(rng);  // avoid zero coeff
        coeffs[i] = double(c);
    }
    return coeffs;
}

/// Random roots in [-spread, spread] with fixed seed
inline std::vector<double> random_roots(int n, double spread = 5.0, uint64_t seed = 20240101) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(-spread, spread);
    std::vector<double> roots(n);
    for (int i = 0; i < n; ++i)
        roots[i] = dist(rng);
    std::sort(roots.begin(), roots.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });
    return roots;
}

// ════════════════════════════════════════════════════════════════
//  Evaluation
// ════════════════════════════════════════════════════════════════

/// Standard Horner: P(t) = sum coeffs[i] * t^{n-i}
template <typename T>
T peval(const std::vector<T>& coeffs, T t) {
    T r = T(0.0);
    for (const auto& c : coeffs) r = r * t + c;
    return r;
}

/// Cauchy bound: |root| <= 1 + max(|a_0|, ..., |a_{n-1}|)
template <typename T>
T cauchy_bound(const std::vector<T>& coeffs) {
    T mx = T(0.0);
    for (size_t i = 1; i < coeffs.size(); ++i) {
        T abs_c = myabs(coeffs[i]);
        if (abs_c > mx) mx = abs_c;
    }
    return T(1.0) + mx;
}

/// Vieta bound: at least one root satisfies |r| <= |a_0|^{1/n}
/// Proof: \prod_{i=1}^n |r_i| = |a_0|.
/// If all |r_i| > |a_0|^{1/n}, the product would exceed |a_0| — contradiction.
/// This bound is typically much tighter than Cauchy.
/// Specialised overloads for double vs mpf_class (nth root differs).
inline double vieta_bound(const std::vector<double>& coeffs) {
    if (coeffs.size() <= 2) return myabs(coeffs.back());
    double a0 = myabs(coeffs.back());
    int n = int(coeffs.size()) - 1;
    return std::pow(a0, 1.0 / n);
}

inline mpf_class vieta_bound(const std::vector<mpf_class>& coeffs) {
    if (coeffs.size() <= 2) return myabs(coeffs.back());
    // Convert to double for rough bound estimate — sufficient as interval endpoint
    double a0_d = myabs(coeffs.back()).get_d();
    int n = int(coeffs.size()) - 1;
    return mpf_class(std::pow(a0_d, 1.0 / n));
}

// ════════════════════════════════════════════════════════════════
//  Discriminant Phi(x) = x * f_{n-1}(x) - a_0
//  Recurrence: f_0 = 1, f_k = a_{n-k} - x * f_{k-1}
//  Returns (Phi(x), quotient [1, f_1, ..., f_{n-1}])
//  IMPORTANT: return value INCLUDES f_0=1 so it can be used
//  directly as the deflated quotient coefficients.
// ════════════════════════════════════════════════════════════════

template <typename T>
std::pair<T, std::vector<T>> phi_eval(const std::vector<T>& coeffs, T x) {
    int n = int(coeffs.size()) - 1;
    std::vector<T> f;           // will be [1, f_1, ..., f_{n-1}]
    f.reserve(n);
    f.push_back(T(1.0));          // f_0 = 1
    for (int k = 1; k < n; ++k) {
        // f_k = a_{n-k} - x * f_{k-1}
        // coeffs[k] = a_{n-k} because coeffs[0]=a_n, coeffs[1]=a_{n-1}, ...
        f.push_back(coeffs[k] - x * f.back());
    }
    // Phi(x) = x * f_{n-1} - a_0
    T phi;
    if (n > 0) {
        phi = x * f.back() - coeffs[n];  // f.back() = f_{n-1}
    } else {
        phi = x * T(1.0) - coeffs[n];
    }
    return {phi, std::move(f)};
}

// ════════════════════════════════════════════════════════════════
//  Standard Horner evaluation: P(x) = a_0 + x*(a_1 + x*(...))
// ════════════════════════════════════════════════════════════════

template <typename T>
T horner_eval(const std::vector<T>& coeffs, T x) {
    int n = int(coeffs.size()) - 1;
    T val = coeffs[0];  // a_n = 1 (monic)
    for (int i = 1; i <= n; ++i) {
        val = val * x + coeffs[i];
    }
    return val;
}

// ════════════════════════════════════════════════════════════════
//  Asymptotic sign for large |x| (avoids overflow)
// ════════════════════════════════════════════════════════════════

/// Sign of P(t) for large |t|, assuming monic polynomial
/// P(t) ~ t^n, so sign = sign(t)^n for large t
inline int asymptotic_sign_P(int n, double t) {
    if (t > 0) return 1;
    return (n % 2 == 0) ? 1 : -1;
}

/// Sign of Phi(x) for large |x|
/// Phi(x) = -P(-x), so sign = -sign((-x)^n)
inline int asymptotic_sign_Phi(int n, double x) {
    if (x > 0) {
        // sign = -sign((-x)^n) = -((-1)^n) = -1 if n even, +1 if n odd
        return (n % 2 == 0) ? -1 : 1;
    } else {
        // sign = -sign((-x)^n) = -(1) = -1 (since -x > 0 for x < 0)
        return -1;
    }
}

// ════════════════════════════════════════════════════════════════
//  Synthetic division (used in HB-ED)
//  Returns quotient coeffs [q_0=1, q_1, ..., q_{n-1}]
// ════════════════════════════════════════════════════════════════

template <typename T>
std::vector<T> synthetic_div(const std::vector<T>& coeffs, T root) {
    int n = int(coeffs.size()) - 1;
    std::vector<T> quot(n, T(0.0));
    T b = coeffs[0];  // a_n
    for (int i = 0; i < n; ++i) {
        quot[i] = b;
        b = b * root + coeffs[i + 1];
    }
    return quot;
}

// ════════════════════════════════════════════════════════════════
//  Residual: max |P(r)| over found roots
// ════════════════════════════════════════════════════════════════

template <typename T>
T max_residual(const std::vector<T>& orig_coeffs, const std::vector<T>& roots) {
    T mx = T(0.0);
    for (const auto& r : roots) {
        T val = myabs(peval(orig_coeffs, r));
        if (val > mx) mx = val;
    }
    return mx;
}

// ════════════════════════════════════════════════════════════════
//  GMP precision helpers
// ════════════════════════════════════════════════════════════════

/// Set GMP default precision (bits) and return previous value
inline mp_bitcnt_t set_gmp_precision(mp_bitcnt_t bits) {
    mp_bitcnt_t old = mpf_get_default_prec();
    mpf_set_default_prec(bits);
    return old;
}

/// Convert vector<double> to vector<mpf_class>
inline std::vector<mpf_class> to_mpf_vec(const std::vector<double>& src) {
    std::vector<mpf_class> dst(src.size());
    for (size_t i = 0; i < src.size(); ++i)
        dst[i] = mpf_class(src[i]);
    return dst;
}

/// Convert vector<mpf_class> to vector<double> (for output)
inline std::vector<double> to_double_vec(const std::vector<mpf_class>& src) {
    std::vector<double> dst(src.size());
    for (size_t i = 0; i < src.size(); ++i)
        dst[i] = src[i].get_d();
    return dst;
}

} // namespace hbd
