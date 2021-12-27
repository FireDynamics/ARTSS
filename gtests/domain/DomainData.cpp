/// \file      DomainData.cpp
/// \brief
/// \date      Dec 16, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include <gtest/gtest.h>
#include "src/utility/settings/Settings.h"
#include "src/domain/DomainData.h"

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
    // config
    Settings::Settings settings;
    settings.sset("domain_parameters/nx", "10");
    settings.sset("domain_parameters/ny", "11");
    settings.sset("domain_parameters/nz", "12");
    settings.sset("domain_parameters/enable_computational_domain", "No");
    settings.sset("domain_parameters/X1", "0");
    settings.sset("domain_parameters/Y1", "-1");
    settings.sset("domain_parameters/Z1", "-0.3");
    settings.sset("domain_parameters/X2", "1");
    settings.sset("domain_parameters/Y2", "1.2");
    settings.sset("domain_parameters/Z2", "0.3");
    settings.sset("solver/description", "non-solver");
    DomainData::init(settings);

    auto domain_data = DomainData::getInstance();

    // size_t
    EXPECT_EQ(10, domain_data->get_nx());
    EXPECT_EQ(11, domain_data->get_ny());
    EXPECT_EQ(12, domain_data->get_nz());

    EXPECT_EQ(12, domain_data->get_Nx());
    EXPECT_EQ(13, domain_data->get_Ny());
    EXPECT_EQ(14, domain_data->get_Nz());

    EXPECT_DOUBLE_EQ( 0.0, domain_data->get_X1());
    EXPECT_DOUBLE_EQ(-1.0, domain_data->get_Y1());
    EXPECT_DOUBLE_EQ(-0.3, domain_data->get_Z1());
    EXPECT_DOUBLE_EQ(1.0, domain_data->get_X2());
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_Y2());
    EXPECT_DOUBLE_EQ(0.3, domain_data->get_Z2());

    EXPECT_DOUBLE_EQ( 0.0, domain_data->get_x1());
    EXPECT_DOUBLE_EQ(-1.0, domain_data->get_y1());
    EXPECT_DOUBLE_EQ(-0.3, domain_data->get_z1());
    EXPECT_DOUBLE_EQ(1.0, domain_data->get_x2());
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_y2());
    EXPECT_DOUBLE_EQ(0.3, domain_data->get_z2());

    EXPECT_DOUBLE_EQ(1.0, domain_data->get_Lx());
    EXPECT_DOUBLE_EQ(2.2, domain_data->get_Ly());
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_Lz());

    EXPECT_DOUBLE_EQ(1.0, domain_data->get_lx());
    EXPECT_DOUBLE_EQ(2.2, domain_data->get_ly());
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_lz());

    EXPECT_DOUBLE_EQ(0.1, domain_data->get_dx());
    EXPECT_DOUBLE_EQ(0.2, domain_data->get_dy());
    EXPECT_DOUBLE_EQ(0.05, domain_data->get_dz());

    EXPECT_EQ(0, domain_data->get_levels());
    EXPECT_EQ(2184, domain_data->get_size());
}

// test with computational domain but same values
TEST_F(DomainDataTest, constructor2Test) {
    // config
    Settings::Settings settings;
    settings.sset("domain_parameters/nx", "10");
    settings.sset("domain_parameters/ny", "11");
    settings.sset("domain_parameters/nz", "12");
    settings.sset("domain_parameters/enable_computational_domain", "Yes");
    settings.sset("domain_parameters/x1", "0");
    settings.sset("domain_parameters/y1", "-1");
    settings.sset("domain_parameters/z1", "-0.3");
    settings.sset("domain_parameters/x2", "1");
    settings.sset("domain_parameters/y2", "1.2");
    settings.sset("domain_parameters/z2", "0.3");
    settings.sset("domain_parameters/X1", "0");
    settings.sset("domain_parameters/Y1", "-1");
    settings.sset("domain_parameters/Z1", "-0.3");
    settings.sset("domain_parameters/X2", "1");
    settings.sset("domain_parameters/Y2", "1.2");
    settings.sset("domain_parameters/Z2", "0.3");
    settings.sset("solver/description", "non-solver");
    DomainData::init(settings);

    auto domain_data = DomainData::getInstance();

    // size_t
    EXPECT_EQ(10, domain_data->get_nx());
    EXPECT_EQ(11, domain_data->get_ny());
    EXPECT_EQ(12, domain_data->get_nz());

    EXPECT_EQ(12, domain_data->get_Nx());
    EXPECT_EQ(13, domain_data->get_Ny());
    EXPECT_EQ(14, domain_data->get_Nz());

    EXPECT_DOUBLE_EQ( 0.0, domain_data->get_X1());
    EXPECT_DOUBLE_EQ(-1.0, domain_data->get_Y1());
    EXPECT_DOUBLE_EQ(-0.3, domain_data->get_Z1());
    EXPECT_DOUBLE_EQ(1.0, domain_data->get_X2());
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_Y2());
    EXPECT_DOUBLE_EQ(0.3, domain_data->get_Z2());

    EXPECT_DOUBLE_EQ( 0.0, domain_data->get_x1());
    EXPECT_DOUBLE_EQ(-1.0, domain_data->get_y1());
    EXPECT_DOUBLE_EQ(-0.3, domain_data->get_z1());
    EXPECT_DOUBLE_EQ(1.0, domain_data->get_x2());
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_y2());
    EXPECT_DOUBLE_EQ(0.3, domain_data->get_z2());

    EXPECT_DOUBLE_EQ(1.0, domain_data->get_Lx());
    EXPECT_DOUBLE_EQ(2.2, domain_data->get_Ly());
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_Lz());

    EXPECT_DOUBLE_EQ(1.0, domain_data->get_lx());
    EXPECT_DOUBLE_EQ(2.2, domain_data->get_ly());
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_lz());

    EXPECT_DOUBLE_EQ(0.1, domain_data->get_dx());
    EXPECT_DOUBLE_EQ(0.2, domain_data->get_dy());
    EXPECT_DOUBLE_EQ(0.05, domain_data->get_dz());

    EXPECT_EQ(0, domain_data->get_levels());
    EXPECT_EQ(2184, domain_data->get_size());
}

// test with computational domain values different from the physical domain
TEST_F(DomainDataTest, constructor3Test) {
    // config
    Settings::Settings settings;
    settings.sset("domain_parameters/nx", "20");
    settings.sset("domain_parameters/ny", "24");
    settings.sset("domain_parameters/nz", "26");
    settings.sset("domain_parameters/enable_computational_domain", "Yes");
    settings.sset("domain_parameters/x1", "0.25");
    settings.sset("domain_parameters/y1", "-0.6");
    settings.sset("domain_parameters/z1", "-0.03");
    settings.sset("domain_parameters/x2", "0.75");
    settings.sset("domain_parameters/y2", "0.6");
    settings.sset("domain_parameters/z2", "0.1");
    settings.sset("domain_parameters/X1", "-0.5");
    settings.sset("domain_parameters/Y1", "-2.2");
    settings.sset("domain_parameters/Z1", "-0.9");
    settings.sset("domain_parameters/X2", "1.5");
    settings.sset("domain_parameters/Y2", "2.2");
    settings.sset("domain_parameters/Z2", "0.3");
    settings.sset("solver/description", "non-solver");
    DomainData::init(settings);

    auto domain_data = DomainData::getInstance();

    // size_t
    EXPECT_EQ(20, domain_data->get_nx());
    EXPECT_EQ(24, domain_data->get_ny());
    EXPECT_EQ(26, domain_data->get_nz());

    EXPECT_EQ(22, domain_data->get_Nx());
    EXPECT_EQ(26, domain_data->get_Ny());
    EXPECT_EQ(28, domain_data->get_Nz());

    EXPECT_DOUBLE_EQ(-0.5, domain_data->get_X1());
    EXPECT_DOUBLE_EQ(-2.2, domain_data->get_Y1());
    EXPECT_DOUBLE_EQ(-0.9, domain_data->get_Z1());
    EXPECT_DOUBLE_EQ(1.5, domain_data->get_X2());
    EXPECT_DOUBLE_EQ(2.2, domain_data->get_Y2());
    EXPECT_DOUBLE_EQ(0.3, domain_data->get_Z2());

    EXPECT_DOUBLE_EQ( 0.25, domain_data->get_x1());
    EXPECT_DOUBLE_EQ(-0.6 , domain_data->get_y1());
    EXPECT_DOUBLE_EQ(-0.03, domain_data->get_z1());
    EXPECT_DOUBLE_EQ(0.75, domain_data->get_x2());
    EXPECT_DOUBLE_EQ(0.6 , domain_data->get_y2());
    EXPECT_DOUBLE_EQ(0.1 , domain_data->get_z2());

    EXPECT_DOUBLE_EQ(2.0, domain_data->get_Lx());
    EXPECT_DOUBLE_EQ(4.4, domain_data->get_Ly());
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_Lz());

    EXPECT_DOUBLE_EQ(0.5 , domain_data->get_lx());
    EXPECT_DOUBLE_EQ(1.2 , domain_data->get_ly());
    EXPECT_DOUBLE_EQ(0.13, domain_data->get_lz());

    EXPECT_DOUBLE_EQ(0.025, domain_data->get_dx());
    EXPECT_DOUBLE_EQ(0.05, domain_data->get_dy());
    EXPECT_DOUBLE_EQ(0.005, domain_data->get_dz());

    EXPECT_EQ(0, domain_data->get_levels());
    EXPECT_EQ(16016, domain_data->get_size());
}

TEST_F(DomainDataTest, goodCaseMultigridTest) {
    // config
    Settings::Settings settings;
    settings.sset("domain_parameters/nx", "24");
    settings.sset("domain_parameters/ny", "4");
    settings.sset("domain_parameters/nz", "1");
    settings.sset("domain_parameters/enable_computational_domain", "No");
    settings.sset("domain_parameters/X1", "0");
    settings.sset("domain_parameters/Y1", "-1");
    settings.sset("domain_parameters/Z1", "-0.3");
    settings.sset("domain_parameters/X2", "1.2");
    settings.sset("domain_parameters/Y2", "-0.2");
    settings.sset("domain_parameters/Z2", "0.3");
    settings.sset("solver/description", "NSSolver");
    settings.sset("solver/pressure/n_level", "2");
    DomainData::init(settings);

    auto domain_data = DomainData::getInstance();

    EXPECT_EQ(24, domain_data->get_nx(0));
    EXPECT_EQ(4, domain_data->get_ny(0));
    EXPECT_EQ(1, domain_data->get_nz(0));

    EXPECT_EQ(12, domain_data->get_nx(1));
    EXPECT_EQ(2, domain_data->get_ny(1));
    EXPECT_EQ(1, domain_data->get_nz(1));

    EXPECT_EQ(6, domain_data->get_nx(2));
    EXPECT_EQ(1, domain_data->get_ny(2));
    EXPECT_EQ(1, domain_data->get_nz(2));

    EXPECT_DOUBLE_EQ(1.2, domain_data->get_lx());
    EXPECT_DOUBLE_EQ(0.8, domain_data->get_ly());
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_lz());

    EXPECT_DOUBLE_EQ(0.05, domain_data->get_dx(0));
    EXPECT_DOUBLE_EQ(0.2, domain_data->get_dy(0));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_dz(0));

    EXPECT_DOUBLE_EQ(0.1, domain_data->get_dx(1));
    EXPECT_DOUBLE_EQ(0.4, domain_data->get_dy(1));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_dz(1));

    EXPECT_DOUBLE_EQ(0.2, domain_data->get_dx(2));
    EXPECT_DOUBLE_EQ(0.8, domain_data->get_dy(2));
    EXPECT_DOUBLE_EQ(0.6, domain_data->get_dz(2));

    EXPECT_EQ(2, domain_data->get_levels());
}

TEST_F(DomainDataTest, badCaseMultigridTest) {
    // config
    Settings::Settings settings;
    settings.sset("domain_parameters/nx", "10");
    settings.sset("domain_parameters/ny", "11");
    settings.sset("domain_parameters/nz", "12");
    settings.sset("domain_parameters/enable_computational_domain", "No");
    settings.sset("domain_parameters/X1", "0");
    settings.sset("domain_parameters/Y1", "-1");
    settings.sset("domain_parameters/Z1", "-0.3");
    settings.sset("domain_parameters/X2", "1");
    settings.sset("domain_parameters/Y2", "1.2");
    settings.sset("domain_parameters/Z2", "0.3");
    settings.sset("solver/description", "NSSolver");
    settings.sset("solver/pressure/n_level", "2");
    DomainData::init(settings);

    auto domain_data = DomainData::getInstance();

    EXPECT_DOUBLE_EQ(0.1, domain_data->get_dx());
    EXPECT_DOUBLE_EQ(0.2, domain_data->get_dy());
    EXPECT_DOUBLE_EQ(0.05, domain_data->get_dz());

    EXPECT_EQ(2, domain_data->get_levels());
    //TODO currently only warnings
}

// bad case: computational domain > physical domain
TEST_F(DomainDataTest, badCaseComputationlDomainTest) {
    // config
    Settings::Settings settings;
    settings.sset("domain_parameters/nx", "20");
    settings.sset("domain_parameters/ny", "22");
    settings.sset("domain_parameters/nz", "24");
    settings.sset("domain_parameters/enable_computational_domain", "Yes");
    settings.sset("domain_parameters/X1", "0.25");
    settings.sset("domain_parameters/Y1", "-0.6");
    settings.sset("domain_parameters/Z1", "-0.03");
    settings.sset("domain_parameters/X2", "0.75");
    settings.sset("domain_parameters/Y2", "0.6");
    settings.sset("domain_parameters/Z2", "0.1");
    settings.sset("domain_parameters/x1", "-0.5");
    settings.sset("domain_parameters/y1", "-2.2");
    settings.sset("domain_parameters/z1", "-0.9");
    settings.sset("domain_parameters/x2", "1.5");
    settings.sset("domain_parameters/y2", "2.2");
    settings.sset("domain_parameters/z2", "0.3");
    settings.sset("solver/description", "non-solver");
    DomainData::init(settings);

    auto domain_data = DomainData::getInstance();

    // size_t
    EXPECT_EQ(20, domain_data->get_nx());
    EXPECT_EQ(22, domain_data->get_ny());
    EXPECT_EQ(24, domain_data->get_nz());

    EXPECT_EQ(22, domain_data->get_Nx());
    EXPECT_EQ(24, domain_data->get_Ny());
    EXPECT_EQ(26, domain_data->get_Nz());

    EXPECT_DOUBLE_EQ(-0.5, domain_data->get_x1());
    EXPECT_DOUBLE_EQ(-2.2, domain_data->get_y1());
    EXPECT_DOUBLE_EQ(-0.9, domain_data->get_z1());
    EXPECT_DOUBLE_EQ(1.5, domain_data->get_x2());
    EXPECT_DOUBLE_EQ(2.2, domain_data->get_y2());
    EXPECT_DOUBLE_EQ(0.3, domain_data->get_z2());

    EXPECT_DOUBLE_EQ( 0.25, domain_data->get_X1());
    EXPECT_DOUBLE_EQ(-0.6 , domain_data->get_Y1());
    EXPECT_DOUBLE_EQ(-0.03, domain_data->get_Z1());
    EXPECT_DOUBLE_EQ(0.75, domain_data->get_X2());
    EXPECT_DOUBLE_EQ(0.6 , domain_data->get_Y2());
    EXPECT_DOUBLE_EQ(0.1 , domain_data->get_Z2());

    EXPECT_DOUBLE_EQ(2.0, domain_data->get_lx());
    EXPECT_DOUBLE_EQ(4.4, domain_data->get_ly());
    EXPECT_DOUBLE_EQ(1.2, domain_data->get_lz());

    EXPECT_DOUBLE_EQ(0.5 , domain_data->get_Lx());
    EXPECT_DOUBLE_EQ(1.2 , domain_data->get_Ly());
    EXPECT_DOUBLE_EQ(0.13, domain_data->get_Lz());

    EXPECT_DOUBLE_EQ(0.1, domain_data->get_dx());
    EXPECT_DOUBLE_EQ(0.2, domain_data->get_dy());
    EXPECT_DOUBLE_EQ(0.05, domain_data->get_dz());

    EXPECT_EQ(0, domain_data->get_levels());

    //TODO something has to happen, computational_domain > physical domain is not allowed
}
