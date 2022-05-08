/// \file       Environment.cpp
/// \brief      
/// \date       Dec 10, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.


#ifndef BENCH_ENVIRONMENT_H_
#define BENCH_ENVIRONMENT_H_

#include <benchmark/benchmark.h>
#include "src/utility/Utility.h"

static int setup_done = 0;

// Override this to define how to set up the environment.
static void DoSetup(const benchmark::State& state) {
    if (setup_done)
        return;

    Utility::create_logger("info", "gtest.log");
    setup_done = 1;
}

#endif /* BENCH_FIELD_FIELD_H_ */
