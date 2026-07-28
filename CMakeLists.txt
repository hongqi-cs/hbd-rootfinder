cmake_minimum_required(VERSION 3.14)
project(HBD_RootFinder VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")

# Find GMP (manual — GMP does not ship CMake config)
find_path(GMP_INCLUDE_DIR gmp.h PATHS /usr/include /usr/local/include)
find_library(GMP_LIBRARY gmp PATHS /usr/lib/x86_64-linux-gnu /usr/lib /usr/local/lib)
find_library(GMPXX_LIBRARY gmpxx PATHS /usr/lib/x86_64-linux-gnu /usr/lib /usr/local/lib)

if(NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
    message(FATAL_ERROR "GMP not found — install libgmp-dev")
endif()

message(STATUS "GMP include: ${GMP_INCLUDE_DIR}")
message(STATUS "GMP libs:    ${GMP_LIBRARY} ; ${GMPXX_LIBRARY}")

# ─── Core library (header-only) ───────────────────────────────
add_library(hbd_core INTERFACE)
target_include_directories(hbd_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src ${GMP_INCLUDE_DIR})
target_link_libraries(hbd_core INTERFACE ${GMP_LIBRARY} ${GMPXX_LIBRARY})

# ─── Experiments ──────────────────────────────────────────────
function(add_hbd_experiment NAME)
    add_executable(${NAME} experiments/${NAME}.cpp)
    target_link_libraries(${NAME} PRIVATE hbd_core)
endfunction()

add_hbd_experiment(exp8_baseline)
add_hbd_experiment(exp4_speed)
add_hbd_experiment(exp5_wilkinson)
add_hbd_experiment(exp6_multiple)
add_hbd_experiment(exp7_memory)
add_hbd_experiment(exp9_vieta)
add_hbd_experiment(run_all)
