/// \file       Coordinate.cpp
/// \brief      
/// \date       Dec 14, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>
#include <src/domain/Coordinate.h>

TEST(CoordinateTest, constructorVal) {
    Coordinate<size_t> coord;
    for (size_t i = 0; i < number_of_axes; i++) {
        EXPECT_EQ(coord[i], 0);
    }
}

TEST(CoordinateTest, constructorVal2) {
    Coordinate<size_t> coord(1, 2, 3);
    for (size_t i = 0; i < number_of_axes; i++) {
        EXPECT_EQ(coord[i], i + 1);
    }
}

TEST(CoordinateTest, fdiff) {
    Coordinate<real> coord1(0, 0, 0);
    Coordinate<real> coord2(1, 2, 1);

    EXPECT_EQ(fdiff(coord2, coord1), coord2);
}

TEST(CoordinateTest, diff) {
    Coordinate<size_t> coord1(0, 0, 0);
    Coordinate<size_t> coord2(1, 2, 1);

    EXPECT_EQ(diff(coord1, coord2), coord2);
}
