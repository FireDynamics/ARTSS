/// \file      DomainData.cpp
/// \brief
/// \date      Dec 16, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include <gtest/gtest.h>
#include "src/utility/settings/Settings.h"
#include "src/domain/DomainData.h"
#include "src/interfaces/ISolver.h"

class DomainDataTest : public testing::Test {
    void SetUp() override {
        DomainData::reset();
    }

    void TearDown() override {
        DomainData::reset();
    }
};

// test without computational domain
TEST_F(DomainDataTest, constructorTest) {
    Settings::domain_parameters domain_parameters{};
    domain_parameters.enable_computational_domain = false;
    domain_parameters.number_of_inner_cells.set_coordinate(10, 11, 12);
    domain_parameters.start_coords_PD.set_coordinate(0, -1, -0.3);
    domain_parameters.end_coords_PD.set_coordinate(1, 1.2, 0.3);
    domain_parameters.start_coords_CD.copy(domain_parameters.start_coords_PD);
    domain_parameters.end_coords_CD.copy(domain_parameters.end_coords_PD);
    // config
    DomainData::init({}, domain_parameters, 0);

    auto domain_data = DomainData::getInstance();

    // size_t
    EXPECT_EQ(10, domain_data->get_number_of_inner_cells(CoordinateAxis::X));
    EXPECT_EQ(11, domain_data->get_number_of_inner_cells(CoordinateAxis::Y));
    EXPECT_EQ(12, domain_data->get_number_of_inner_cells(CoordinateAxis::Z));

    EXPECT_EQ(12, domain_data->get_number_of_cells(CoordinateAxis::X));
    EXPECT_EQ(13, domain_data->get_number_of_cells(CoordinateAxis::Y));
    EXPECT_EQ(14, domain_data->get_number_of_cells(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(0.0, domain_data->get_start_coord_PD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(-1.0, domain_data->get_start_coord_PD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(-0.3, domain_data->get_start_coord_PD(CoordinateAxis::Z));
    EXPECT_DOUBLE_EQ(1.0, domain_data->get_end_coord_PD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_end_coord_PD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.3, domain_data->get_end_coord_PD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(0.0, domain_data->get_start_coord_CD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(-1.0, domain_data->get_start_coord_CD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(-0.3, domain_data->get_start_coord_CD(CoordinateAxis::Z));
    EXPECT_DOUBLE_EQ(1.0, domain_data->get_end_coord_CD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_end_coord_CD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.3, domain_data->get_end_coord_CD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(1.0, domain_data->get_length_PD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(2.2, domain_data->get_length_PD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_length_PD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(1.0, domain_data->get_length_CD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(2.2, domain_data->get_length_CD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_length_CD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(0.1, domain_data->get_spacing(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(0.2, domain_data->get_spacing(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.05, domain_data->get_spacing(CoordinateAxis::Z));

    EXPECT_EQ(0, domain_data->get_levels());
    EXPECT_EQ(2184, domain_data->get_size());
}

// test with computational domain but same values
TEST_F(DomainDataTest, constructor2Test) {
    Settings::domain_parameters domain_parameters{};
    domain_parameters.enable_computational_domain = true;
    domain_parameters.number_of_inner_cells.set_coordinate(10, 11, 12);
    domain_parameters.start_coords_PD.set_coordinate(0, -1, -0.3);
    domain_parameters.end_coords_PD.set_coordinate(1, 1.2, 0.3);
    domain_parameters.start_coords_CD.copy(domain_parameters.start_coords_PD);
    domain_parameters.end_coords_CD.copy(domain_parameters.end_coords_PD);
    // config
    DomainData::init({}, domain_parameters, 0);

    auto domain_data = DomainData::getInstance();

    // size_t
    EXPECT_EQ(10, domain_data->get_number_of_inner_cells(CoordinateAxis::X));
    EXPECT_EQ(11, domain_data->get_number_of_inner_cells(CoordinateAxis::Y));
    EXPECT_EQ(12, domain_data->get_number_of_inner_cells(CoordinateAxis::Z));

    EXPECT_EQ(12, domain_data->get_number_of_cells(CoordinateAxis::X));
    EXPECT_EQ(13, domain_data->get_number_of_cells(CoordinateAxis::Y));
    EXPECT_EQ(14, domain_data->get_number_of_cells(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(0.0, domain_data->get_start_coord_PD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(-1.0, domain_data->get_start_coord_PD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(-0.3, domain_data->get_start_coord_PD(CoordinateAxis::Z));
    EXPECT_DOUBLE_EQ(1.0, domain_data->get_end_coord_PD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_end_coord_PD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.3, domain_data->get_end_coord_PD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(0.0, domain_data->get_start_coord_CD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(-1.0, domain_data->get_start_coord_CD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(-0.3, domain_data->get_start_coord_CD(CoordinateAxis::Z));
    EXPECT_DOUBLE_EQ(1.0, domain_data->get_end_coord_CD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_end_coord_CD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.3, domain_data->get_end_coord_CD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(1.0, domain_data->get_length_PD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(2.2, domain_data->get_length_PD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_length_PD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(1.0, domain_data->get_length_CD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(2.2, domain_data->get_length_CD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_length_CD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(0.1, domain_data->get_spacing(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(0.2, domain_data->get_spacing(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.05, domain_data->get_spacing(CoordinateAxis::Z));

    EXPECT_EQ(0, domain_data->get_levels());
    EXPECT_EQ(2184, domain_data->get_size());
}

// test with computational domain values different from the physical domain
TEST_F(DomainDataTest, constructor3Test) {
    // config
    Settings::domain_parameters domain_parameters{};
    domain_parameters.enable_computational_domain = true;
    domain_parameters.number_of_inner_cells.set_coordinate(20, 24, 26);
    domain_parameters.start_coords_CD.set_coordinate(0.25, -0.6, -0.03);
    domain_parameters.end_coords_CD.set_coordinate(0.75, 0.6, 0.1);
    domain_parameters.start_coords_PD.set_coordinate(-0.5, -2.2, -0.9);
    domain_parameters.end_coords_PD.set_coordinate(1.5, 2.2, 0.3);
    DomainData::init({}, domain_parameters, 0);

    auto domain_data = DomainData::getInstance();

    // size_t
    EXPECT_EQ(20, domain_data->get_number_of_inner_cells(CoordinateAxis::X));
    EXPECT_EQ(24, domain_data->get_number_of_inner_cells(CoordinateAxis::Y));
    EXPECT_EQ(26, domain_data->get_number_of_inner_cells(CoordinateAxis::Z));

    EXPECT_EQ(22, domain_data->get_number_of_cells(CoordinateAxis::X));
    EXPECT_EQ(26, domain_data->get_number_of_cells(CoordinateAxis::Y));
    EXPECT_EQ(28, domain_data->get_number_of_cells(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(-0.5, domain_data->get_start_coord_PD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(-2.2, domain_data->get_start_coord_PD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(-0.9, domain_data->get_start_coord_PD(CoordinateAxis::Z));
    EXPECT_DOUBLE_EQ(1.5, domain_data->get_end_coord_PD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(2.2, domain_data->get_end_coord_PD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.3, domain_data->get_end_coord_PD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(0.25, domain_data->get_start_coord_CD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(-0.6, domain_data->get_start_coord_CD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(-0.03, domain_data->get_start_coord_CD(CoordinateAxis::Z));
    EXPECT_DOUBLE_EQ(0.75, domain_data->get_end_coord_CD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_end_coord_CD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.1, domain_data->get_end_coord_CD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(2.0, domain_data->get_length_PD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(4.4, domain_data->get_length_PD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_length_PD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(0.5, domain_data->get_length_CD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_length_CD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.13, domain_data->get_length_CD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(0.025, domain_data->get_spacing(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(0.05, domain_data->get_spacing(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.005, domain_data->get_spacing(CoordinateAxis::Z));

    EXPECT_EQ(0, domain_data->get_levels());
    EXPECT_EQ(16016, domain_data->get_size());
}

TEST_F(DomainDataTest, goodCaseMultigridTest) {
    Settings::domain_parameters domain_parameters{};
    domain_parameters.enable_computational_domain = false;
    domain_parameters.number_of_inner_cells.set_coordinate(24, 4, 1);
    domain_parameters.start_coords_PD.set_coordinate(0, -1, -0.3);
    domain_parameters.end_coords_PD.set_coordinate(1.2, -0.2, 0.3);
    domain_parameters.start_coords_CD.copy(domain_parameters.start_coords_PD);
    domain_parameters.end_coords_CD.copy(domain_parameters.end_coords_PD);
    DomainData::init({}, domain_parameters, 2);
    auto domain_data = DomainData::getInstance();

    EXPECT_EQ(24, domain_data->get_number_of_inner_cells(CoordinateAxis::X, 0));
    EXPECT_EQ(4 , domain_data->get_number_of_inner_cells(CoordinateAxis::Y, 0));
    EXPECT_EQ(1 , domain_data->get_number_of_inner_cells(CoordinateAxis::Z, 0));

    EXPECT_EQ(12, domain_data->get_number_of_inner_cells(CoordinateAxis::X, 1));
    EXPECT_EQ(2 , domain_data->get_number_of_inner_cells(CoordinateAxis::Y, 1));
    EXPECT_EQ(1 , domain_data->get_number_of_inner_cells(CoordinateAxis::Z, 1));

    EXPECT_EQ(6, domain_data->get_number_of_inner_cells(CoordinateAxis::X, 2));
    EXPECT_EQ(1, domain_data->get_number_of_inner_cells(CoordinateAxis::Y, 2));
    EXPECT_EQ(1, domain_data->get_number_of_inner_cells(CoordinateAxis::Z, 2));

    EXPECT_EQ(26, domain_data->get_number_of_cells(CoordinateAxis::X, 0));
    EXPECT_EQ(6 , domain_data->get_number_of_cells(CoordinateAxis::Y, 0));
    EXPECT_EQ(3 , domain_data->get_number_of_cells(CoordinateAxis::Z, 0));

    EXPECT_EQ(14, domain_data->get_number_of_cells(CoordinateAxis::X, 1));
    EXPECT_EQ(4 , domain_data->get_number_of_cells(CoordinateAxis::Y, 1));
    EXPECT_EQ(3 , domain_data->get_number_of_cells(CoordinateAxis::Z, 1));

    EXPECT_EQ(8, domain_data->get_number_of_cells(CoordinateAxis::X, 2));
    EXPECT_EQ(3, domain_data->get_number_of_cells(CoordinateAxis::Y, 2));
    EXPECT_EQ(3, domain_data->get_number_of_cells(CoordinateAxis::Z, 2));

    EXPECT_DOUBLE_EQ(1.2, domain_data->get_length_CD(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(0.8, domain_data->get_length_CD(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_length_CD(CoordinateAxis::Z));

    EXPECT_DOUBLE_EQ(0.05, domain_data->get_spacing(CoordinateAxis::X, 0));
    EXPECT_DOUBLE_EQ(0.2, domain_data->get_spacing(CoordinateAxis::Y, 0));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_spacing(CoordinateAxis::Z, 0));

    EXPECT_DOUBLE_EQ(0.1, domain_data->get_spacing(CoordinateAxis::X, 1));
    EXPECT_DOUBLE_EQ(0.4, domain_data->get_spacing(CoordinateAxis::Y, 1));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_spacing(CoordinateAxis::Z, 1));

    EXPECT_DOUBLE_EQ(0.2, domain_data->get_spacing(CoordinateAxis::X, 2));
    EXPECT_DOUBLE_EQ(0.8, domain_data->get_spacing(CoordinateAxis::Y, 2));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_spacing(CoordinateAxis::Z, 2));

    EXPECT_EQ(2, domain_data->get_levels());
}

TEST_F(DomainDataTest, badCaseMultigridTest) {
    Settings::domain_parameters domain_parameters{};
    domain_parameters.enable_computational_domain = false;
    domain_parameters.number_of_inner_cells.set_coordinate(10, 11, 12);
    domain_parameters.start_coords_PD.set_coordinate(0, -1, -0.3);
    domain_parameters.end_coords_PD.set_coordinate(1, 1.2, 0.3);
    domain_parameters.start_coords_CD.copy(domain_parameters.start_coords_PD);
    domain_parameters.end_coords_CD.copy(domain_parameters.end_coords_PD);
    DomainData::init({}, domain_parameters, 2);
    auto domain_data = DomainData::getInstance();

    EXPECT_DOUBLE_EQ(0.1, domain_data->get_spacing(CoordinateAxis::X));
    EXPECT_DOUBLE_EQ(0.2, domain_data->get_spacing(CoordinateAxis::Y));
    EXPECT_DOUBLE_EQ(0.05, domain_data->get_spacing(CoordinateAxis::Z));

    EXPECT_EQ(2, domain_data->get_levels());
    //TODO currently only warnings
}
