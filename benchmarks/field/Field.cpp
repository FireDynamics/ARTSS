#include <benchmark/benchmark.h>

#include "../Environment.h"
#include "src/field/Field.h"

#include <math.h>


static void BM_AddScalar(benchmark::State &state) {
    // Utility::create_logger("info", "bench.log");
    size_t size = state.range(0);

    for (auto _ : state) {
        Field a(UNKNOWN_FIELD, 0.0, 0, size);
        a += M_E;
    }
}

static void BM_AddFields(benchmark::State &state) {
    // Utility::create_logger("info", "bench.log");
    size_t size = state.range(0);

    for (auto _ : state) {
        Field a(UNKNOWN_FIELD, 0.0, 0, size);
        Field b(UNKNOWN_FIELD, 0.0, 0, size);
        a += b;
    }
}

static void BM_MulScalar(benchmark::State &state) {
    // Utility::create_logger("info", "bench.log");
    size_t size = state.range(0);

    for (auto _ : state) {
        Field a(UNKNOWN_FIELD, 0.0, 0, size);
        a *= M_E;
    }
}

static void BM_MulFields(benchmark::State &state) {
    // Utility::create_logger("info", "bench.log");
    size_t size = state.range(0);

    for (auto _ : state) {
        Field a(UNKNOWN_FIELD, 0.0, 0, size);
        Field b(UNKNOWN_FIELD, 0.0, 0, size);
        a *= b;
    }
}

BENCHMARK(BM_AddScalar)->ThreadRange(1, 8)->Range(8, 8<<20)->Setup(DoSetup);
BENCHMARK(BM_AddFields)->ThreadRange(1, 8)->Range(8, 8<<20)->Setup(DoSetup);
BENCHMARK(BM_MulScalar)->ThreadRange(1, 8)->Range(8, 8<<20)->Setup(DoSetup);
BENCHMARK(BM_MulFields)->ThreadRange(1, 8)->Range(8, 8<<20)->Setup(DoSetup);

BENCHMARK_MAIN();
