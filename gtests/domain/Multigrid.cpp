/// \file      Multigrid.cpp
/// \brief
/// \date      Dec 31, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include <gtest/gtest.h>
#include "src/utility/settings/Settings.h"
#include "src/domain/Obstacle.h"
#include "src/domain/Surface.h"
#include "src/domain/Multigrid.h"

class MultigridTest : public testing::Test {
    void SetUp() override {
        Settings::domain_parameters domain_params{};
        domain_params.start_coords_PD.set_coordinate(0, 0, 0);
        domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
        domain_params.end_coords_PD.set_coordinate(5, 10, 15);
        domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
        domain_params.enable_computational_domain = false;
        domain_params.number_of_inner_cells.set_coordinate(8, 12, 4);
        DomainData::init({}, domain_params, 2);
    }

    void TearDown() override {
        DomainData::reset();
    }
};

TEST_F(MultigridTest, sizes) {
    std::vector<Surface> surfaces;
    std::vector<BoundaryDataController> bdc_surfaces;
    std::vector<Obstacle> obstacles;
    std::vector<BoundaryDataController> bdc_obstacles;
    BoundaryDataController bdc_domain({});
    const size_t level = 2;

    Multigrid multigrid(surfaces, bdc_surfaces, obstacles, bdc_obstacles, bdc_domain, level);

    EXPECT_EQ(840, multigrid.get_slice_size_domain_cells(0));
    EXPECT_EQ(384, multigrid.get_slice_size_domain_inner_cells_level_joined(0));

    EXPECT_EQ(192, multigrid.get_slice_size_domain_cells(1));
    EXPECT_EQ(48, multigrid.get_slice_size_domain_inner_cells_level_joined(1));

    EXPECT_EQ(60, multigrid.get_slice_size_domain_cells(2));
    EXPECT_EQ(6, multigrid.get_slice_size_domain_inner_cells_level_joined(2));

    EXPECT_EQ(1092, multigrid.get_size_domain_cells());
    EXPECT_EQ(438, multigrid.get_size_domain_inner_cells_level_joined());

    EXPECT_EQ(0, multigrid.get_slice_size_obstacle_cells(0));
    EXPECT_EQ(0, multigrid.get_slice_size_obstacle_cells(1));
    EXPECT_EQ(0, multigrid.get_slice_size_obstacle_cells(2));
    EXPECT_EQ(0, multigrid.get_size_obstacle_cells());
}