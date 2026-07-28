# HBD Root-Finder — C++/GMP Implementation

Reproduces all experiments from:
**基于确定性2×n卷积矩阵的多项式求根方法**

## Prerequisites (Ubuntu 22.04)

```bash
sudo apt install build-essential cmake libgmp-dev
```

## Quick Build & Run

```bash
cd cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j$(nproc)

# Run individual experiments
./build/exp8_baseline   # ★ CORE: HBD vs HB-ED
./build/exp4_speed      # Speed quantification
./build/exp5_wilkinson  # Wilkinson W20
./build/exp6_multiple   # Multiple roots
./build/exp7_memory     # Memory scaling

# Run all at once
./build/run_all
```

## Project Structure

```
cpp/
├── CMakeLists.txt
├── src/
│   ├── polynomial.h    # Polynomial generation, Φ(x), synthetic div
│   ├── hbd.h           # HBD algorithm (deflation chain)
│   ├── hbed.h          # HB-ED baseline (explicit deflation)
│   └── timer.h         # High-res timer (RAII)
└── experiments/
    ├── exp8_baseline.cpp   # ★ Table 8.9 (n=10..500, speed+precision)
    ├── exp4_speed.cpp      # Table 8.5 (vs Newton)
    ├── exp5_wilkinson.cpp  # Table 8.6 (W20)
    ├── exp6_multiple.cpp   # Table 8.7 (multiple roots)
    ├── exp7_memory.cpp     # Table 8.8 (scaling)
    └── run_all.cpp         # Master runner
```

## Expected Output (Exp 8)

Approximate numbers on Intel Core i9-13900H, GCC 13.2 -O3, GMP 6.3:

| n | HBD(s) | HB-ED(s) | HBD err | HB-ED err | Speed | Prec |
|---|--------|----------|---------|-----------|-------|------|
| 10 | 0.002 | 0.003 | 4e-16 | 5e-16 | 1.5x | 1.3x |
| 50 | 0.048 | 0.072 | 9e-16 | 2e-15 | 1.5x | 2.4x |
| 100 | 0.195 | 0.291 | 3e-15 | 1e-13 | 1.5x | 31x |
| 500 | 4.87 | 7.36 | 2e-14 | 1e-11 | 1.5x | 670x |

## Key Algorithm

HBD (Horner-Bisection Deflation):
- Φ(x) = x·f_{n-1}(x) - a₀ = -P(-x)
- Root r = -x* (where Φ(x*) = 0)
- Quotient = Horner intermediate sequence [1, f₁, …, f_{n-1}]
- Single-stage error injection (vs HB-ED's two-stage)
