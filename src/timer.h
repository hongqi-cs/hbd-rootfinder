#pragma once
#include <chrono>
#include <string>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

namespace hbd {

/// High-resolution wall-clock timer
class Timer {
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = Clock::time_point;
    TimePoint t0;
public:
    Timer() : t0(Clock::now()) {}
    void reset() { t0 = Clock::now(); }
    /// Elapsed seconds
    double secs() const {
        auto d = Clock::now() - t0;
        return std::chrono::duration<double>(d).count();
    }
    /// Elapsed milliseconds
    double msecs() const { return secs() * 1000.0; }
};

/// Timing helper: run fn() n_times, return avg/max/min seconds
template <typename F>
struct TimingStats {
    double avg, max, min;
};
template <typename F>
TimingStats<double> time_it(F fn, int n_runs = 10) {
    std::vector<double> ts(n_runs);
    for (int i = 0; i < n_runs; ++i) {
        Timer t;
        fn();
        ts[i] = t.secs();
    }
    double avg = std::accumulate(ts.begin(), ts.end(), 0.0) / n_runs;
    auto [mn, mx] = std::minmax_element(ts.begin(), ts.end());
    return {avg, *mx, *mn};
}

/// RAII scoped timer: prints label + elapsed when destroyed
class ScopedTimer {
    std::string label;
    Timer timer;
public:
    explicit ScopedTimer(std::string lbl) : label(std::move(lbl)) {}
    ~ScopedTimer() {
        std::cout << "  [" << label << "] " << timer.secs() << " s\n";
    }
};

} // namespace hbd
