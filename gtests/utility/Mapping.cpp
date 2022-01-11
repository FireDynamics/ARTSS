/// \file       Mapping.cpp
/// \brief      
/// \date       Dec 14, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>
#include "src/domain/Coordinate.h"

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

TEST(MappingTest, getBoundaryPatchName) {
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

TEST(MappingTest, getAxisName) {
    EXPECT_EQ("X", Mapping::get_axis_name(X));
    EXPECT_EQ("Y", Mapping::get_axis_name(Y));
    EXPECT_EQ("Z", Mapping::get_axis_name(Z));
}

TEST(MappingTest, matchField) {
    FieldType field_type_rho = Mapping::match_field("rho");
    ASSERT_EQ(field_type_rho, FieldType::RHO);
    FieldType field_type_u = Mapping::match_field("u");
    ASSERT_EQ(field_type_u, FieldType::U);
    FieldType field_type_v = Mapping::match_field("v");
    ASSERT_EQ(field_type_v, FieldType::V);
    FieldType field_type_w = Mapping::match_field("w");
    ASSERT_EQ(field_type_w, FieldType::W);
    FieldType field_type_p = Mapping::match_field("p");
    ASSERT_EQ(field_type_p, FieldType::P);
    FieldType field_type_T = Mapping::match_field("T");
    ASSERT_EQ(field_type_T, FieldType::T);
    FieldType field_type_unknown = Mapping::match_field("a");
    ASSERT_EQ(field_type_unknown, FieldType::UNKNOWN_FIELD);
}

TEST(MappingTest, matchBoundaryCondition) {
    BoundaryCondition bc_neumann = Mapping::match_boundary_condition("neumann");
    ASSERT_EQ(bc_neumann, BoundaryCondition::NEUMANN);
    BoundaryCondition bc_dirichlet = Mapping::match_boundary_condition("dirichlet");
    ASSERT_EQ(bc_dirichlet, BoundaryCondition::DIRICHLET);
    BoundaryCondition bc_periodic = Mapping::match_boundary_condition("periodic");
    ASSERT_EQ(bc_periodic, BoundaryCondition::PERIODIC);
    BoundaryCondition bc_unknown = Mapping::match_boundary_condition("a");
    ASSERT_EQ(bc_unknown, BoundaryCondition::UNKNOWN_CONDITION);
}

TEST(MappingTest, getFieldTypeName) {
    std::string rho = Mapping::get_field_type_name(FieldType::RHO);
    ASSERT_EQ("rho", rho);
    std::string u = Mapping::get_field_type_name(FieldType::U);
    ASSERT_EQ("u", u);
    std::string v = Mapping::get_field_type_name(FieldType::V);
    ASSERT_EQ("v", v);
    std::string w = Mapping::get_field_type_name(FieldType::W);
    ASSERT_EQ("w", w);
    std::string p = Mapping::get_field_type_name(FieldType::P);
    ASSERT_EQ("p", p);
    std::string T = Mapping::get_field_type_name(FieldType::T);
    ASSERT_STRCASEEQ("T", T.c_str());

//    ASSERT_EXIT((Mapping::get_field_type_name(FieldType::UNKNOWN_FIELD),exit(0)),::testing::KilledBySignal(SIGEV_SIGNAL),".*");
}

TEST(MappingTest, toAxis) {
    EXPECT_EQ(CoordinateAxis::X, Mapping::to_axis(Patch::LEFT));
    EXPECT_EQ(CoordinateAxis::X, Mapping::to_axis(Patch::RIGHT));
    EXPECT_EQ(CoordinateAxis::Y, Mapping::to_axis(Patch::BOTTOM));
    EXPECT_EQ(CoordinateAxis::Y, Mapping::to_axis(Patch::TOP));
    EXPECT_EQ(CoordinateAxis::Z, Mapping::to_axis(Patch::FRONT));
    EXPECT_EQ(CoordinateAxis::Z, Mapping::to_axis(Patch::BACK));
}

TEST(MappingTest, toPatch) {
    EXPECT_EQ(Patch::LEFT, Mapping::to_patch(CoordinateAxis::X, true));
    EXPECT_EQ(Patch::RIGHT, Mapping::to_patch(CoordinateAxis::X, false));
    EXPECT_EQ(Patch::BOTTOM, Mapping::to_patch(CoordinateAxis::Y, true));
    EXPECT_EQ(Patch::TOP, Mapping::to_patch(CoordinateAxis::Y, false));
    EXPECT_EQ(Patch::FRONT, Mapping::to_patch(CoordinateAxis::Z, true));
    EXPECT_EQ(Patch::BACK, Mapping::to_patch(CoordinateAxis::Z, false));
}
