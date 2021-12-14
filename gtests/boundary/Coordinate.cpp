/// \file       Coordinate.cpp
/// \brief      
/// \date       Dec 14, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>
#include <src/boundary/Coordinate.h>
#include "src/boundary/PatchObject.h"

TEST(CoordinateTest, constructorVal) {
    Coordinate<size_t> coord;
    for (size_t i = 0; i < number_of_axes; i++) {
        EXPECT_EQ(coord[i], 0);
    }
}

TEST(CoordinateTest, constructorVal2) {
    Coordinate<size_t> coord(1,2,3);
    for (size_t i = 0; i < number_of_axes; i++) {
        EXPECT_EQ(coord[i], i+1);
    }
}

TEST(CoordinateTest, matchCoordinate) {
    EXPECT_EQ(X, Axis::match_axis("X"));
    EXPECT_EQ(Y, Axis::match_axis("Y"));
    EXPECT_EQ(Z, Axis::match_axis("Z"));
    EXPECT_EQ(X, Axis::match_axis("x"));
    EXPECT_EQ(Y, Axis::match_axis("y"));
    EXPECT_EQ(Z, Axis::match_axis("z"));
}

TEST(CoordinateTest, getName) {
    EXPECT_EQ("X", Axis::get_axis_name(X));
    EXPECT_EQ("Y", Axis::get_axis_name(Y));
    EXPECT_EQ("Z", Axis::get_axis_name(Z));
}

