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

    // domain cells
    EXPECT_EQ(840, multigrid.get_slice_size_domain_cells(0));
    EXPECT_EQ(192, multigrid.get_slice_size_domain_cells(1));
    EXPECT_EQ(60, multigrid.get_slice_size_domain_cells(2));
    EXPECT_EQ(1092, multigrid.get_size_domain_cells());

    // domain inner cells
    EXPECT_EQ(384, multigrid.get_slice_size_domain_inner_cells_level_joined(0));
    EXPECT_EQ(48, multigrid.get_slice_size_domain_inner_cells_level_joined(1));
    EXPECT_EQ(6, multigrid.get_slice_size_domain_inner_cells_level_joined(2));
    EXPECT_EQ(438, multigrid.get_size_domain_inner_cells_level_joined());

    // obstacle cells
    EXPECT_EQ(0, multigrid.get_slice_size_obstacle_cells(0));
    EXPECT_EQ(0, multigrid.get_slice_size_obstacle_cells(1));
    EXPECT_EQ(0, multigrid.get_slice_size_obstacle_cells(2));
    EXPECT_EQ(0, multigrid.get_size_obstacle_cells());
}

TEST_F(MultigridTest, indices) {
    std::vector<Surface> surfaces;
    std::vector<BoundaryDataController> bdc_surfaces;
    std::vector<Obstacle> obstacles;
    std::vector<BoundaryDataController> bdc_obstacles;
    BoundaryDataController bdc_domain({});
    const size_t level = 2;
    Multigrid multigrid(surfaces, bdc_surfaces, obstacles, bdc_obstacles, bdc_domain, level);

    // obstacle cells
    EXPECT_EQ(0, multigrid.get_start_index_obstacle_cells_level_joined(0));
    EXPECT_EQ(0, multigrid.get_end_index_obstacle_cells_level_joined(0));

    EXPECT_EQ(0, multigrid.get_start_index_obstacle_cells_level_joined(1));
    EXPECT_EQ(0, multigrid.get_end_index_obstacle_cells_level_joined(1));

    EXPECT_EQ(0, multigrid.get_start_index_obstacle_cells_level_joined(2));
    EXPECT_EQ(0, multigrid.get_end_index_obstacle_cells_level_joined(2));

    // domain cells
    EXPECT_EQ(0, multigrid.get_start_index_domain_cells_level_joined(0));
    EXPECT_EQ(839, multigrid.get_end_index_domain_cells_level_joined(0));

    EXPECT_EQ(840, multigrid.get_start_index_domain_cells_level_joined(1));
    EXPECT_EQ(1031, multigrid.get_end_index_domain_cells_level_joined(1));

    EXPECT_EQ(1032, multigrid.get_start_index_domain_cells_level_joined(2));
    EXPECT_EQ(1091, multigrid.get_end_index_domain_cells_level_joined(2));

    // domain inner cells
    EXPECT_EQ(0, multigrid.get_start_index_domain_inner_cells_level_joined(0));
    EXPECT_EQ(383, multigrid.get_end_index_domain_inner_cells_level_joined(0));

    EXPECT_EQ(384, multigrid.get_start_index_domain_inner_cells_level_joined(1));
    EXPECT_EQ(431, multigrid.get_end_index_domain_inner_cells_level_joined(1));

    EXPECT_EQ(432, multigrid.get_start_index_domain_inner_cells_level_joined(2));
    EXPECT_EQ(437, multigrid.get_end_index_domain_inner_cells_level_joined(2));
}

TEST_F(MultigridTest, Obstacle_even_size) {
    DomainData::reset();

    Settings::domain_parameters domain_params{};
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(8, 8, 8);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(64, 64, 64);
    DomainData::init({}, domain_params, 2);


    std::vector<Surface> surfaces;
    std::vector<BoundaryDataController> bdc_surfaces;
    Coordinate<size_t> o_start(25, 25, 25);
    Coordinate<size_t> o_end(40, 40, 40);
    std::vector<Obstacle> obstacles;
    obstacles.reserve(1);
    obstacles.emplace_back(o_start, o_end, 0, "test_obstacle");
    std::vector<BoundaryDataController> bdc_obstacles;
    BoundaryDataController bdc_domain({});
    const size_t level = 2;
    Multigrid multigrid(surfaces, bdc_surfaces, obstacles, bdc_obstacles, bdc_domain, level);

    // obstacle cells

    // level 0
    EXPECT_EQ(0, multigrid.get_start_index_obstacle_cells_level_joined(0));
    EXPECT_EQ(4095, multigrid.get_end_index_obstacle_cells_level_joined(0));

    EXPECT_EQ(25+25*66+25*66*66, multigrid.get_obstacle_cells_level_joined()[0]);
    EXPECT_EQ(40+40*66+40*66*66, multigrid.get_obstacle_cells_level_joined()[4095]);

    // level 1
    EXPECT_EQ(4096, multigrid.get_start_index_obstacle_cells_level_joined(1));
    EXPECT_EQ(4607, multigrid.get_end_index_obstacle_cells_level_joined(1));

    EXPECT_EQ(13+13*34+13*34*34, multigrid.get_obstacle_cells_level_joined()[4096]);
    EXPECT_EQ(20+20*34+20*34*34, multigrid.get_obstacle_cells_level_joined()[4607]);

    // level 2
    EXPECT_EQ(4608, multigrid.get_start_index_obstacle_cells_level_joined(2));
    EXPECT_EQ(4671, multigrid.get_end_index_obstacle_cells_level_joined(2));

    EXPECT_EQ(7+7*18+7*18*18, multigrid.get_obstacle_cells_level_joined()[4608]);
    EXPECT_EQ(10+10*18+10*18*18, multigrid.get_obstacle_cells_level_joined()[4671]);
}

TEST_F(MultigridTest, Obstacle_odd_size) {
    DomainData::reset();

    Settings::domain_parameters domain_params{};
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(8, 8, 8);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(64, 64, 64);
    DomainData::init({}, domain_params, 2);


    std::vector<Surface> surfaces;
    std::vector<BoundaryDataController> bdc_surfaces;
    Coordinate<size_t> o_start(25, 25, 25);
    Coordinate<size_t> o_end(41, 40, 40);
    std::vector<Obstacle> obstacles;
    obstacles.reserve(1);
    obstacles.emplace_back(o_start, o_end, 0, "test_obstacle");
    std::vector<BoundaryDataController> bdc_obstacles;
    BoundaryDataController bdc_domain({});
    const size_t level = 2;
    Multigrid multigrid(surfaces, bdc_surfaces, obstacles, bdc_obstacles, bdc_domain, level);

    // obstacle cells

    // level 0
    EXPECT_EQ(0, multigrid.get_start_index_obstacle_cells_level_joined(0));
    EXPECT_EQ(4351, multigrid.get_end_index_obstacle_cells_level_joined(0));

    EXPECT_EQ(25+25*66+25*66*66, multigrid.get_obstacle_cells_level_joined()[0]);
    EXPECT_EQ(41+40*66+40*66*66, multigrid.get_obstacle_cells_level_joined()[4095]);

    // level 1
    EXPECT_EQ(4352, multigrid.get_start_index_obstacle_cells_level_joined(1));
    EXPECT_EQ(4927, multigrid.get_end_index_obstacle_cells_level_joined(1));

    EXPECT_EQ(13+13*34+13*34*34, multigrid.get_obstacle_cells_level_joined()[4096]);
    EXPECT_EQ(21+20*34+20*34*34, multigrid.get_obstacle_cells_level_joined()[4607]);

    // level 2
    EXPECT_EQ(4928, multigrid.get_start_index_obstacle_cells_level_joined(2));
    EXPECT_EQ(5007, multigrid.get_end_index_obstacle_cells_level_joined(2));

    EXPECT_EQ(7+7*18+7*18*18, multigrid.get_obstacle_cells_level_joined()[4608]);
    EXPECT_EQ(11+10*18+10*18*18, multigrid.get_obstacle_cells_level_joined()[4671]);
}
