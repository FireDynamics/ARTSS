/// \file       BoundaryData.cpp
/// \brief      
/// \date       Sep 02, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include <gtest/gtest.h>
#include <src/boundary/BoundaryData.h>

class BoundaryDataTest : public testing::Test {
};

TEST_F(BoundaryDataTest, match_field) {
    FieldType field_type_rho = BoundaryData::match_field("rho");
    ASSERT_EQ(field_type_rho, FieldType::RHO);
    FieldType field_type_u = BoundaryData::match_field("u");
    ASSERT_EQ(field_type_u, FieldType::U);
    FieldType field_type_v = BoundaryData::match_field("v");
    ASSERT_EQ(field_type_v, FieldType::V);
    FieldType field_type_w = BoundaryData::match_field("w");
    ASSERT_EQ(field_type_w, FieldType::W);
    FieldType field_type_p = BoundaryData::match_field("p");
    ASSERT_EQ(field_type_p, FieldType::P);
    FieldType field_type_T = BoundaryData::match_field("T");
    ASSERT_EQ(field_type_T, FieldType::T);
    FieldType field_type_unknown = BoundaryData::match_field("a");
    ASSERT_EQ(field_type_unknown, FieldType::UNKNOWN_FIELD);
}

TEST_F(BoundaryDataTest, match_patch) {
    Patch patch_front = BoundaryData::match_patch("front");
    ASSERT_EQ(patch_front, Patch::FRONT);
    Patch patch_back = BoundaryData::match_patch("back");
    ASSERT_EQ(patch_back, Patch::BACK);
    Patch patch_bottom = BoundaryData::match_patch("bottom");
    ASSERT_EQ(patch_bottom, Patch::BOTTOM);
    Patch patch_top = BoundaryData::match_patch("top");
    ASSERT_EQ(patch_top, Patch::TOP);
    Patch patch_left = BoundaryData::match_patch("left");
    ASSERT_EQ(patch_left, Patch::LEFT);
    Patch patch_right = BoundaryData::match_patch("right");
    ASSERT_EQ(patch_right, Patch::RIGHT);
    Patch patch_unknown = BoundaryData::match_patch("a");
    ASSERT_EQ(patch_unknown, Patch::UNKNOWN_PATCH);
}

TEST_F(BoundaryDataTest, match_boundary_condition) {
    BoundaryCondition bc_neumann = BoundaryData::match_boundary_condition("neumann");
    ASSERT_EQ(bc_neumann, BoundaryCondition::NEUMANN);
    BoundaryCondition bc_dirichlet = BoundaryData::match_boundary_condition("dirichlet");
    ASSERT_EQ(bc_dirichlet, BoundaryCondition::DIRICHLET);
    BoundaryCondition bc_periodic = BoundaryData::match_boundary_condition("periodic");
    ASSERT_EQ(bc_periodic, BoundaryCondition::PERIODIC);
    BoundaryCondition bc_unknown = BoundaryData::match_boundary_condition("a");
    ASSERT_EQ(bc_unknown, BoundaryCondition::UNKNOWN_CONDITION);
}