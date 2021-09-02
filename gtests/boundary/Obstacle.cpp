/// \file       Obstacle.cpp
/// \brief      
/// \date       Sep 02, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>
#include <src/boundary/Obstacle.h>

class ObstacleTest : public testing::Test {
};

TEST_F(ObstacleTest, constructor_val) {
    real x1 = 0;
    real x2 = 1;
    real y1 = 0;
    real y2 = 1;
    real z1 = 0;
    real z2 = 1;
    Obstacle *obstacle = new Obstacle(x1, x2, y1, y2, z1, z2, "test obstacle");
}
