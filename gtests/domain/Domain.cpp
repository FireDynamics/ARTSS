/// \file      Domain.cpp
/// \brief
/// \date      Dec 31, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include <gtest/gtest.h>
#include "src/utility/settings/Settings.h"
#include "src/domain/Domain.h"
#include "src/domain/Obstacle.h"
#include "src/domain/Surface.h"

TEST(DomainTest, sizes) {
    Settings::domain_parameters domain_params {};
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({}, domain_params, 0);

    size_t size_obstacle_list = 0;
    auto obstacle_list = new size_t[size_obstacle_list];
    auto surface_list = new size_t*[number_of_patches];
    PatchObject size_surface_list;
    size_t level = 0;
    Domain domain(obstacle_list, size_obstacle_list, surface_list, size_surface_list, level);

    PatchObject *sizes = domain.get_size_boundary_list();
    PatchObject boundary_sizes(*sizes);
    EXPECT_EQ(12*22, boundary_sizes[Patch::LEFT]);
    EXPECT_EQ(12*22, boundary_sizes[Patch::RIGHT]);
    EXPECT_EQ(7*22, boundary_sizes[Patch::FRONT]);
    EXPECT_EQ(7*22, boundary_sizes[Patch::BACK]);
    EXPECT_EQ(7*12, boundary_sizes[Patch::BOTTOM]);
    EXPECT_EQ(7*12, boundary_sizes[Patch::TOP]);
    EXPECT_EQ(1004, boundary_sizes.get_sum());
    EXPECT_EQ(1848, domain.get_size_domain_list());
    EXPECT_EQ(1000, domain.get_size_inner_list());
}

TEST(DomainTest, sizesWithObstacle) {
    Settings::domain_parameters domain_params {};
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({}, domain_params, 0);

    size_t size_obstacle_list = 9;
    size_t obstacle_list[] = {162, 163, 164, 169, 170, 171, 176, 177, 178};
    auto surface_list = new size_t*[number_of_patches];
    PatchObject size_surface_list;
    size_t level = 0;
    Domain domain(obstacle_list, size_obstacle_list, surface_list, size_surface_list, level);

    PatchObject *sizes = domain.get_size_boundary_list();
    PatchObject boundary_sizes(*sizes);
    EXPECT_EQ(12*22, boundary_sizes[Patch::LEFT]);
    EXPECT_EQ(12*22, boundary_sizes[Patch::RIGHT]);
    EXPECT_EQ(7*22, boundary_sizes[Patch::FRONT]);
    EXPECT_EQ(7*22, boundary_sizes[Patch::BACK]);
    EXPECT_EQ(7*12, boundary_sizes[Patch::BOTTOM]);
    EXPECT_EQ(7*12, boundary_sizes[Patch::TOP]);
    EXPECT_EQ(1004, boundary_sizes.get_sum());
    EXPECT_EQ(1839, domain.get_size_domain_list());
    EXPECT_EQ(991, domain.get_size_inner_list());

    delete[] surface_list;
}

TEST(DomainTest, sizesWithSurfaces) {
    Settings::domain_parameters domain_params {};
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({}, domain_params, 0);

    size_t size_obstacle_list = 0;
    size_t obstacle_list[] = {};
    PatchObject size_surface_list;
    size_surface_list[Patch::LEFT] = 3;
    size_surface_list[Patch::FRONT] = 1;
    auto surface_list = new size_t*[number_of_patches];
    for (size_t p = 0; p < number_of_patches; p++) {
        surface_list[p] = new size_t[size_surface_list[p]];
    }
    surface_list[Patch::LEFT][0] = 322;
    surface_list[Patch::LEFT][1] = 476;
    surface_list[Patch::LEFT][2] = 630;
    surface_list[Patch::FRONT][0] = 33;
    size_t level = 0;
    Domain domain(obstacle_list, size_obstacle_list, surface_list, size_surface_list, level);

    PatchObject *sizes = domain.get_size_boundary_list();
    PatchObject boundary_sizes(*sizes);
    EXPECT_EQ(261, boundary_sizes[Patch::LEFT]);
    EXPECT_EQ(264, boundary_sizes[Patch::RIGHT]);
    EXPECT_EQ(153, boundary_sizes[Patch::FRONT]);
    EXPECT_EQ(154, boundary_sizes[Patch::BACK]);
    EXPECT_EQ(84, boundary_sizes[Patch::BOTTOM]);
    EXPECT_EQ(84, boundary_sizes[Patch::TOP]);
    EXPECT_EQ(1000, boundary_sizes.get_sum());
    EXPECT_EQ(1844, domain.get_size_domain_list());
    EXPECT_EQ(1000, domain.get_size_inner_list());
    delete[] surface_list;
}

TEST(DomainTest, sizesWithSurfaceAndObstacle) {
    Settings::domain_parameters domain_params {};
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({}, domain_params, 0);

    size_t size_obstacle_list = 9;
    size_t obstacle_list[] = {162, 163, 164, 169, 170, 171, 176, 177, 178};
    PatchObject size_surface_list;
    size_surface_list[Patch::LEFT] = 3;
    size_surface_list[Patch::FRONT] = 1;
    auto surface_list = new size_t*[number_of_patches];
    for (size_t p = 0; p < number_of_patches; p++) {
        surface_list[p] = new size_t[size_surface_list[p]];
    }
    surface_list[Patch::LEFT][0] = 322;
    surface_list[Patch::LEFT][1] = 476;
    surface_list[Patch::LEFT][2] = 630;
    surface_list[Patch::FRONT][0] = 33;
    size_t level = 0;
    Domain domain(obstacle_list, size_obstacle_list, surface_list, size_surface_list, level);

    PatchObject *sizes = domain.get_size_boundary_list();
    PatchObject boundary_sizes(*sizes);
    EXPECT_EQ(261, boundary_sizes[Patch::LEFT]);
    EXPECT_EQ(264, boundary_sizes[Patch::RIGHT]);
    EXPECT_EQ(153, boundary_sizes[Patch::FRONT]);
    EXPECT_EQ(154, boundary_sizes[Patch::BACK]);
    EXPECT_EQ(84, boundary_sizes[Patch::BOTTOM]);
    EXPECT_EQ(84, boundary_sizes[Patch::TOP]);
    EXPECT_EQ(1000, boundary_sizes.get_sum());
    EXPECT_EQ(1835, domain.get_size_domain_list());
    EXPECT_EQ(991, domain.get_size_inner_list());
    delete[] surface_list;
}

TEST(DomainTest, order) {
    Settings::domain_parameters domain_params{};
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({}, domain_params, 0);

    size_t size_obstacle_list = 0;
    auto obstacle_list = new size_t[size_obstacle_list];
    auto surface_list = new size_t*[number_of_patches];
    PatchObject size_surface_list;
    size_t level = 0;
    Domain domain(obstacle_list, size_obstacle_list, surface_list, size_surface_list, level);

    auto inner_list = domain.get_inner_list();
    for (size_t i = 1; i < domain.get_size_inner_list(); i++) {
        EXPECT_GT(inner_list[i], inner_list[i - 1]);
    }

    auto boundary_lists = domain.get_boundary_list();
    auto boundary_sizes = domain.get_size_boundary_list();
    for (Patch p: all_patches) {
        auto boundary_list = boundary_lists[p];
        size_t boundary_size = (*boundary_sizes)[p];
        for (size_t i = 1; i < boundary_size; i++) {
            EXPECT_GT(boundary_list[i], boundary_list[i - 1]);
        }
    }
}
