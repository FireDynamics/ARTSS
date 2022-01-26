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
    Settings::domain_parameters domain_params{ };
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({ }, domain_params, 0);

    size_t size_obstacle_list = 0;
    std::vector<size_t> obstacle_list(0);
    std::vector<size_t *> surface_list(number_of_patches);
    PatchObject size_surface_list;
    size_t level = 0;
    Domain domain(obstacle_list.data(), obstacle_list.size(), surface_list.data(), size_surface_list, level);

    PatchObject *sizes = domain.get_size_boundary_list();
    PatchObject boundary_sizes(*sizes);
    EXPECT_EQ(12 * 22, boundary_sizes[Patch::LEFT]);
    EXPECT_EQ(12 * 22, boundary_sizes[Patch::RIGHT]);
    EXPECT_EQ(7 * 22, boundary_sizes[Patch::FRONT]);
    EXPECT_EQ(7 * 22, boundary_sizes[Patch::BACK]);
    EXPECT_EQ(7 * 12, boundary_sizes[Patch::BOTTOM]);
    EXPECT_EQ(7 * 12, boundary_sizes[Patch::TOP]);
    EXPECT_EQ(1004, boundary_sizes.get_sum());
    EXPECT_EQ(1848, domain.get_size_domain_list());
    EXPECT_EQ(1000, domain.get_size_inner_list());
}

TEST(DomainTest, sizesWithComputationalDomain) {
    Settings::domain_parameters domain_params{ };
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.set_coordinate(0, 0, 0);
    domain_params.end_coords_CD.set_coordinate(5, 10, 15);
    domain_params.end_coords_PD.set_coordinate(10, 11, 15);
    domain_params.enable_computational_domain = true;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({ }, domain_params, 0);

    size_t size_obstacle_list = 0;
    std::vector<size_t> obstacle_list(0);
    std::vector<size_t *> surface_list(number_of_patches);
    PatchObject size_surface_list;
    size_t level = 0;
    Domain domain(obstacle_list.data(), obstacle_list.size(), surface_list.data(), size_surface_list, level);

    PatchObject *sizes = domain.get_size_boundary_list();
    EXPECT_EQ(12 * 22, (*sizes)[Patch::LEFT]);
    EXPECT_EQ(12 * 22, (*sizes)[Patch::RIGHT]);
    EXPECT_EQ(7 * 22, (*sizes)[Patch::FRONT]);
    EXPECT_EQ(7 * 22, (*sizes)[Patch::BACK]);
    EXPECT_EQ(7 * 12, (*sizes)[Patch::BOTTOM]);
    EXPECT_EQ(7 * 12, (*sizes)[Patch::TOP]);
    EXPECT_EQ(1004, sizes->get_sum());
    EXPECT_EQ(1848, domain.get_size_domain_list());
    EXPECT_EQ(1000, domain.get_size_inner_list());
}

TEST(DomainTest, sizesLevelOne) {
    Settings::domain_parameters domain_params{ };
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(8, 12, 4);
    DomainData::init({ }, domain_params, 1);

    size_t size_obstacle_list = 0;
    std::vector<size_t> obstacle_list(size_obstacle_list);
    std::vector<size_t *> surface_list(number_of_patches);
    PatchObject size_surface_list;
    size_t level = 1;
    Domain domain(obstacle_list.data(), obstacle_list.size(), surface_list.data(), size_surface_list, level);

    PatchObject *sizes = domain.get_size_boundary_list();
    PatchObject boundary_sizes(*sizes);
    EXPECT_EQ(32, boundary_sizes[Patch::LEFT]);
    EXPECT_EQ(32, boundary_sizes[Patch::RIGHT]);
    EXPECT_EQ(24, boundary_sizes[Patch::BOTTOM]);
    EXPECT_EQ(24, boundary_sizes[Patch::TOP]);
    EXPECT_EQ(48, boundary_sizes[Patch::FRONT]);
    EXPECT_EQ(48, boundary_sizes[Patch::BACK]);
    EXPECT_EQ(208, boundary_sizes.get_sum());
    EXPECT_EQ(192, domain.get_size_domain_list());
    EXPECT_EQ(48, domain.get_size_inner_list());
}

TEST(DomainTest, sizesWithObstacle) {
    Settings::domain_parameters domain_params{ };
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({ }, domain_params, 0);

    std::vector<size_t> obstacle_list = {162, 163, 164, 169, 170, 171, 176, 177, 178};
    std::vector<size_t *> surface_list(number_of_patches);
    PatchObject size_surface_list;
    size_t level = 0;
    Domain domain(obstacle_list.data(), obstacle_list.size(), surface_list.data(), size_surface_list, level);

    PatchObject *sizes = domain.get_size_boundary_list();
    PatchObject boundary_sizes(*sizes);
    EXPECT_EQ(12 * 22, boundary_sizes[Patch::LEFT]);
    EXPECT_EQ(12 * 22, boundary_sizes[Patch::RIGHT]);
    EXPECT_EQ(7 * 22, boundary_sizes[Patch::FRONT]);
    EXPECT_EQ(7 * 22, boundary_sizes[Patch::BACK]);
    EXPECT_EQ(7 * 12, boundary_sizes[Patch::BOTTOM]);
    EXPECT_EQ(7 * 12, boundary_sizes[Patch::TOP]);
    EXPECT_EQ(1004, boundary_sizes.get_sum());
    EXPECT_EQ(1839, domain.get_size_domain_list());
    EXPECT_EQ(991, domain.get_size_inner_list());
}

TEST(DomainTest, sizesWithSurfaces) {
    Settings::domain_parameters domain_params{ };
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({ }, domain_params, 0);

    size_t size_obstacle_list = 0;
    std::vector<size_t> obstacle_list(size_obstacle_list);
    PatchObject size_surface_list;
    size_surface_list[Patch::LEFT] = 3;
    size_surface_list[Patch::FRONT] = 1;
    std::vector<size_t *> surface_list(number_of_patches);
    for (size_t p = 0; p < number_of_patches; p++) {
        surface_list[p] = new size_t[size_surface_list[p]];
    }
    surface_list[Patch::LEFT][0] = 322;
    surface_list[Patch::LEFT][1] = 476;
    surface_list[Patch::LEFT][2] = 630;
    surface_list[Patch::FRONT][0] = 33;
    size_t level = 0;
    Domain domain(obstacle_list.data(), obstacle_list.size(), surface_list.data(), size_surface_list, level);

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
}

TEST(DomainTest, sizesWithSurfaceAndObstacle) {
    Settings::domain_parameters domain_params{ };
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({ }, domain_params, 0);

    std::vector<size_t> obstacle_list = {162, 163, 164, 169, 170, 171, 176, 177, 178};
    PatchObject size_surface_list;
    size_surface_list[Patch::LEFT] = 3;
    size_surface_list[Patch::FRONT] = 1;
    std::vector<size_t *> surface_list(number_of_patches);
    for (size_t p = 0; p < number_of_patches; p++) {
        surface_list[p] = new size_t[size_surface_list[p]];
    }
    surface_list[Patch::LEFT][0] = 322;
    surface_list[Patch::LEFT][1] = 476;
    surface_list[Patch::LEFT][2] = 630;
    surface_list[Patch::FRONT][0] = 33;
    size_t level = 0;
    Domain domain(obstacle_list.data(), obstacle_list.size(), surface_list.data(), size_surface_list, level);

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
    for (size_t p = 0; p < number_of_patches; p++) {
        delete[] surface_list[p];
    }
}

TEST(DomainTest, order) {
    Settings::domain_parameters domain_params{ };
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({ }, domain_params, 0);

    size_t size_obstacle_list = 0;
    std::vector<size_t> obstacle_list(size_obstacle_list);
    std::vector<size_t *> surface_list(number_of_patches);
    PatchObject size_surface_list;
    size_t level = 0;
    Domain domain(obstacle_list.data(), obstacle_list.size(), surface_list.data(), size_surface_list, level);

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
}

TEST(DomainTest, indexTestComputationalDomain) {
    Settings::domain_parameters domain_params{ };
    domain_params.start_coords_CD.set_coordinate(0, 0, 0);
    domain_params.start_coords_PD.set_coordinate(0, 0, -1.5);
    domain_params.end_coords_CD.set_coordinate(5, 10, 15);
    domain_params.end_coords_PD.set_coordinate(7, 15, 30);
    domain_params.enable_computational_domain = true;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({ }, domain_params, 0);

    size_t size_obstacle_list = 0;
    std::vector<size_t> obstacle_list(size_obstacle_list);
    std::vector<size_t *> surface_list(number_of_patches);
    PatchObject size_surface_list;
    size_t level = 0;
    Domain domain(obstacle_list.data(), obstacle_list.size(), surface_list.data(), size_surface_list, level);

    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_number_of_cells(X);
    size_t Ny = domain_data->get_number_of_cells(Y);
    size_t start_x = domain_data->get_start_index_CD(CoordinateAxis::X);
    size_t start_y = domain_data->get_start_index_CD(CoordinateAxis::Y);
    size_t start_z = domain_data->get_start_index_CD(CoordinateAxis::Z);
    size_t start_index = IX(start_x, start_y, start_z, Nx, Ny);
    size_t end_x = domain_data->get_end_index_CD(CoordinateAxis::X);
    size_t end_y = domain_data->get_end_index_CD(CoordinateAxis::Y);
    size_t end_z = domain_data->get_end_index_CD(CoordinateAxis::Z);
    size_t end_index = IX(end_x, end_y, end_z, Nx, Ny);

    auto inner_list = domain.get_inner_list();
    EXPECT_EQ(inner_list[0], 586);
    EXPECT_EQ(inner_list[0], start_index);
    EXPECT_EQ(inner_list[domain.get_size_inner_list() - 1], 3353);
    EXPECT_EQ(inner_list[domain.get_size_inner_list() - 1], end_index);

    auto boundary_lists = domain.get_boundary_list();
    auto boundary_sizes = domain.get_size_boundary_list();
    size_t start_boundary = IX(start_x - 1, start_y - 1, start_z - 1, Nx, Ny);
    EXPECT_EQ(boundary_lists[Patch::LEFT][0], start_boundary);
    EXPECT_EQ(boundary_lists[Patch::LEFT][0], 288);
    EXPECT_EQ(boundary_lists[Patch::BOTTOM][0], 288);
    EXPECT_EQ(boundary_lists[Patch::FRONT][0], 288);
    size_t start_right = IX(end_x + 1, start_y - 1, start_z - 1, Nx, Ny);
    EXPECT_EQ(boundary_lists[Patch::RIGHT][0], start_right);
    EXPECT_EQ(boundary_lists[Patch::RIGHT][0], 294);
    size_t start_top = IX(start_x - 1, end_y + 1, start_z - 1, Nx, Ny);
    EXPECT_EQ(boundary_lists[Patch::TOP][0], start_top);
    EXPECT_EQ(boundary_lists[Patch::TOP][0], 477);
    size_t start_back = IX(start_x - 1, start_y - 1, end_z + 1, Nx, Ny);
    EXPECT_EQ(boundary_lists[Patch::BACK][0], start_back);
    EXPECT_EQ(boundary_lists[Patch::BACK][0], 3456);

    size_t end_boundary = IX(end_x + 1, end_y + 1, end_z + 1, Nx, Ny);
    size_t end_left = IX(start_x - 1, end_y + 1, end_z + 1, Nx, Ny);
    EXPECT_EQ(boundary_lists[Patch::LEFT][(*boundary_sizes)[Patch::LEFT] - 1], 3645);
    EXPECT_EQ(boundary_lists[Patch::LEFT][(*boundary_sizes)[Patch::LEFT] - 1], end_left);
    size_t end_bottom = IX(end_x + 1, start_y - 1, end_z + 1, Nx, Ny);
    EXPECT_EQ(boundary_lists[Patch::BOTTOM][(*boundary_sizes)[Patch::BOTTOM] - 1], 3462);
    EXPECT_EQ(boundary_lists[Patch::BOTTOM][(*boundary_sizes)[Patch::BOTTOM] - 1], end_bottom);
    size_t end_front = IX(end_x + 1, end_y + 1, start_z - 1, Nx, Ny);
    EXPECT_EQ(boundary_lists[Patch::FRONT][(*boundary_sizes)[Patch::FRONT] - 1], 483);
    EXPECT_EQ(boundary_lists[Patch::FRONT][(*boundary_sizes)[Patch::FRONT] - 1], end_front);
    EXPECT_EQ(boundary_lists[Patch::RIGHT][(*boundary_sizes)[Patch::RIGHT] - 1], 3651);
    EXPECT_EQ(boundary_lists[Patch::TOP][(*boundary_sizes)[Patch::TOP] - 1], 3651);
    EXPECT_EQ(boundary_lists[Patch::BACK][(*boundary_sizes)[Patch::BACK] - 1], 3651);
    EXPECT_EQ(boundary_lists[Patch::BACK][(*boundary_sizes)[Patch::BACK] - 1], end_boundary);


    auto domain_list = domain.get_domain_list();
    size_t size_domain = domain.get_size_domain_list();
    EXPECT_EQ(288, domain_list[0]);
    EXPECT_EQ(3651, domain_list[size_domain - 1]);
}

TEST(DomainTest, indexTestPhysicalDomain) {
    Settings::domain_parameters domain_params{ };
    domain_params.start_coords_PD.set_coordinate(0, 0, 0);
    domain_params.start_coords_CD.copy(domain_params.start_coords_PD);
    domain_params.end_coords_PD.set_coordinate(5, 10, 15);
    domain_params.end_coords_CD.copy(domain_params.end_coords_PD);
    domain_params.enable_computational_domain = false;
    domain_params.number_of_inner_cells.set_coordinate(5, 20, 10);
    DomainData::init({ }, domain_params, 0);

    size_t size_obstacle_list = 0;
    std::vector<size_t> obstacle_list(size_obstacle_list);
    std::vector<size_t *> surface_list(number_of_patches);
    PatchObject size_surface_list;
    size_t level = 0;
    Domain domain(obstacle_list.data(), obstacle_list.size(), surface_list.data(), size_surface_list, level);

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
}