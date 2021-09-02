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
