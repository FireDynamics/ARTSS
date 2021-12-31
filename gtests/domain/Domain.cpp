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
    delete[] surface_list;
    delete[] obstacle_list;
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

    for (size_t p = 0; p < number_of_patches; p++) {
        delete[] surface_list[p];
    }
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

    auto domain_list = domain.get_domain_list();
    for (size_t i = 1; i < domain.get_size_domain_list(); i++) {
        EXPECT_GT(domain_list[i], domain_list[i - 1]);
    }
    delete[] surface_list;
    delete[] obstacle_list;
}

TEST(DomainTest, indexTestComputationalDomain) {
    Settings::domain_parameters domain_params{};
    domain_params.start_coords_CD.set_coordinate(0, 0, 0);
    domain_params.start_coords_PD.set_coordinate(0, 0, -1.5);
    domain_params.end_coords_CD.set_coordinate(5, 10, 15);
    domain_params.end_coords_PD.set_coordinate(7, 15, 30);
    domain_params.enable_computational_domain = true;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({}, domain_params, 0);

    size_t size_obstacle_list = 0;
    auto obstacle_list = new size_t[size_obstacle_list];
    auto surface_list = new size_t *[number_of_patches];
    PatchObject size_surface_list;
    size_t level = 0;
    Domain domain(obstacle_list, size_obstacle_list, surface_list, size_surface_list, level);

    auto inner_list = domain.get_inner_list();
    EXPECT_EQ(inner_list[0], 316);
    EXPECT_EQ(inner_list[domain.get_size_inner_list() - 1], 1839);

    auto boundary_lists = domain.get_boundary_list();
    auto boundary_sizes = domain.get_size_boundary_list();
    EXPECT_EQ(boundary_lists[Patch::LEFT][0], 154);
    EXPECT_EQ(boundary_lists[Patch::BOTTOM][0], 154);
    EXPECT_EQ(boundary_lists[Patch::FRONT][0], 154);
    EXPECT_EQ(boundary_lists[Patch::RIGHT][0], 160);
    EXPECT_EQ(boundary_lists[Patch::TOP][0], 301);
    EXPECT_EQ(boundary_lists[Patch::BACK][0], 1848);

    EXPECT_EQ(boundary_lists[Patch::LEFT][(*boundary_sizes)[Patch::LEFT] - 1], 1995);
    EXPECT_EQ(boundary_lists[Patch::BOTTOM][(*boundary_sizes)[Patch::BOTTOM] - 1], 1854);
    EXPECT_EQ(boundary_lists[Patch::FRONT][(*boundary_sizes)[Patch::FRONT] - 1], 307);
    EXPECT_EQ(boundary_lists[Patch::RIGHT][(*boundary_sizes)[Patch::RIGHT] - 1], 2001);
    EXPECT_EQ(boundary_lists[Patch::TOP][(*boundary_sizes)[Patch::TOP] - 1], 2001);
    EXPECT_EQ(boundary_lists[Patch::BACK][(*boundary_sizes)[Patch::BACK] - 1], 2001);


    auto domain_list = domain.get_domain_list();
    size_t size_domain = domain.get_size_domain_list();
    EXPECT_EQ(154, domain_list[0]);
    EXPECT_EQ(2001, domain_list[size_domain - 1]);
    delete[] surface_list;
    delete[] obstacle_list;
}

TEST(DomainTest, indexTestPhysicalDomain) {
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
    auto surface_list = new size_t *[number_of_patches];
    PatchObject size_surface_list;
    size_t level = 0;
    Domain domain(obstacle_list, size_obstacle_list, surface_list, size_surface_list, level);

    auto inner_list = domain.get_inner_list();
    EXPECT_EQ(inner_list[0], 162);
    EXPECT_EQ(inner_list[domain.get_size_inner_list() - 1], 1685);

    auto boundary_lists = domain.get_boundary_list();
    auto boundary_sizes = domain.get_size_boundary_list();
    EXPECT_EQ(boundary_lists[Patch::LEFT][0], 0);
    EXPECT_EQ(boundary_lists[Patch::BOTTOM][0], 0);
    EXPECT_EQ(boundary_lists[Patch::FRONT][0], 0);
    EXPECT_EQ(boundary_lists[Patch::RIGHT][0], 6);
    EXPECT_EQ(boundary_lists[Patch::TOP][0], 147);
    EXPECT_EQ(boundary_lists[Patch::BACK][0], 1694);

    EXPECT_EQ(boundary_lists[Patch::LEFT][(*boundary_sizes)[Patch::LEFT] - 1], 1841);
    EXPECT_EQ(boundary_lists[Patch::BOTTOM][(*boundary_sizes)[Patch::BOTTOM] - 1], 1700);
    EXPECT_EQ(boundary_lists[Patch::FRONT][(*boundary_sizes)[Patch::FRONT] - 1], 153);
    EXPECT_EQ(boundary_lists[Patch::RIGHT][(*boundary_sizes)[Patch::RIGHT] - 1], 1847);
    EXPECT_EQ(boundary_lists[Patch::TOP][(*boundary_sizes)[Patch::TOP] - 1], 1847);
    EXPECT_EQ(boundary_lists[Patch::BACK][(*boundary_sizes)[Patch::BACK] - 1], 1847);


    auto domain_list = domain.get_domain_list();
    size_t size_domain = domain.get_size_domain_list();
    EXPECT_EQ(0, domain_list[0]);
    EXPECT_EQ(1847, domain_list[size_domain - 1]);
    delete[] surface_list;
    delete[] obstacle_list;
}
