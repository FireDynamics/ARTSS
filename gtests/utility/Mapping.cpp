/// \file       Mapping.cpp
/// \brief      
/// \date       Dec 14, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>
#include <src/boundary/Coordinate.h>

TEST(MappingTest, matchPatch) {
    Patch patch_front = Mapping::match_patch("front");
    EXPECT_EQ(patch_front, Patch::FRONT);
    Patch patch_back = Mapping::match_patch("back");
    EXPECT_EQ(patch_back, Patch::BACK);
    Patch patch_bottom = Mapping::match_patch("bottom");
    EXPECT_EQ(patch_bottom, Patch::BOTTOM);
    Patch patch_top = Mapping::match_patch("top");
    EXPECT_EQ(patch_top, Patch::TOP);
    Patch patch_left = Mapping::match_patch("left");
    EXPECT_EQ(patch_left, Patch::LEFT);
    Patch patch_right = Mapping::match_patch("right");
    EXPECT_EQ(patch_right, Patch::RIGHT);
    Patch patch_unknown = Mapping::match_patch("a");
    EXPECT_EQ(patch_unknown, Patch::UNKNOWN_PATCH);
}

TEST(MappingTest, get_boundary_patch_name) {
    std::string front = Mapping::get_patch_name(Patch::FRONT);
    EXPECT_EQ("front", front);
    std::string back = Mapping::get_patch_name(Patch::BACK);
    EXPECT_EQ("back", back);
    std::string bottom = Mapping::get_patch_name(Patch::BOTTOM);
    EXPECT_EQ("bottom", bottom);
    std::string top = Mapping::get_patch_name(Patch::TOP);
    EXPECT_EQ("top", top);
    std::string left = Mapping::get_patch_name(Patch::LEFT);
    EXPECT_EQ("left", left);
    std::string right = Mapping::get_patch_name(Patch::RIGHT);
    EXPECT_EQ("right", right);
}

TEST(MappingTest, coordinateAxis) {
    EXPECT_EQ(Mapping::to_patch(CoordinateAxis::X, true), Patch::LEFT);
    EXPECT_EQ(Mapping::to_patch(CoordinateAxis::X, false), Patch::RIGHT);
    EXPECT_EQ(Mapping::to_patch(CoordinateAxis::Y, true), Patch::BOTTOM);
    EXPECT_EQ(Mapping::to_patch(CoordinateAxis::Y, false), Patch::TOP);
    EXPECT_EQ(Mapping::to_patch(CoordinateAxis::Z, true), Patch::FRONT);
    EXPECT_EQ(Mapping::to_patch(CoordinateAxis::Z, false), Patch::BACK);
}

TEST(MappingTest, matchCoordinate) {
    EXPECT_EQ(X, Mapping::match_axis("X"));
    EXPECT_EQ(Y, Mapping::match_axis("Y"));
    EXPECT_EQ(Z, Mapping::match_axis("Z"));
    EXPECT_EQ(X, Mapping::match_axis("x"));
    EXPECT_EQ(Y, Mapping::match_axis("y"));
    EXPECT_EQ(Z, Mapping::match_axis("z"));
}

TEST(MappingTest, getName) {
    EXPECT_EQ("X", Mapping::get_axis_name(X));
    EXPECT_EQ("Y", Mapping::get_axis_name(Y));
    EXPECT_EQ("Z", Mapping::get_axis_name(Z));
}

