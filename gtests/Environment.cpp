/// \file       Environment.cpp
/// \brief      
/// \date       Dec 10, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>
#include "src/utility/Utility.h"

class Environment : public ::testing::Environment {
public:
    ~Environment() override {}

    // Override this to define how to set up the environment.
    void SetUp() override {
        Utility::create_logger("info", "gtest.log");
    }

    // Override this to define how to tear down the environment.
    void TearDown() override {}
};

testing::Environment* const env = testing::AddGlobalTestEnvironment(new Environment);
