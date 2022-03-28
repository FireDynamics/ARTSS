/// \file       SLAdvect.cpp
/// \brief      
/// \date       Mar 22, 2022
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2022> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>
#include <src/advection/SLAdvect.h>
#include "src/domain/DomainData.h"

class SLAdvectTest : public testing::Test {
    void SetUp() override {
        Settings::domain_parameters domain_parameters{ };
        domain_parameters.enable_computational_domain = false;
        domain_parameters.number_of_inner_cells.set_coordinate(10, 11, 12);
        domain_parameters.start_coords_PD.set_coordinate(0, -1, -0.3);
        domain_parameters.end_coords_PD.set_coordinate(1, 1.2, 0.3);
        domain_parameters.start_coords_CD.copy(domain_parameters.start_coords_PD);
        domain_parameters.end_coords_CD.copy(domain_parameters.end_coords_PD);
        DomainData::init({ }, domain_parameters, 0);
    }
};

TEST_F(SLAdvectTest, linear_back_trace2) {
    auto domain_data = DomainData::getInstance();
    real epsilon = 1e-6;
    real trace_back = 1.7;

    std::vector<size_t> solution_i0 = {2, 3, 4, 5, 6, 7, 8, 9, 10, 10};
    std::vector<size_t> solution_i1 = {3, 4, 5, 6, 7, 8, 9, 10, 11, 11};
    real solution_r = 0.7;
    for (size_t i = 1; i <= 10; i++) {
        auto[i0, i1, r] = SLAdvect::calculate_backward_index(CoordinateAxis::X, i, epsilon, trace_back);
        EXPECT_EQ(i0, solution_i0[i - 1]);
        EXPECT_EQ(i1, solution_i1[i - 1]);
        EXPECT_DOUBLE_EQ(r, solution_r);
    }

    std::vector<size_t> solution_j0 = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11};
    std::vector<size_t> solution_j1 = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12};
    real solution_s = 0.7;
    for (size_t j = 1; j <= 11; j++) {
        auto[j0, j1, s] = SLAdvect::calculate_backward_index(CoordinateAxis::Y, j, epsilon, trace_back);
        EXPECT_EQ(j0, solution_j0[j - 1]);
        EXPECT_EQ(j1, solution_j1[j - 1]);
        EXPECT_DOUBLE_EQ(s, solution_s);
    }

    std::vector<size_t> solution_k0 = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12};
    std::vector<size_t> solution_k1 = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13};
    real solution_t = 0.7;
    for (size_t k = 1; k <= 12; k++) {
        auto[k0, k1, t] = SLAdvect::calculate_backward_index(CoordinateAxis::Z, k, epsilon, trace_back);
        EXPECT_EQ(k0, solution_k0[k - 1]);
        EXPECT_EQ(k1, solution_k1[k - 1]);
        EXPECT_DOUBLE_EQ(t, solution_t);
    }
}

TEST_F(SLAdvectTest, linear_back_trace) {
    auto domain_data = DomainData::getInstance();
    real epsilon = 1e-6;
    real trace_back = -1.2;

    std::vector<size_t> solution_i0 = {2, 3, 4, 5, 6, 7, 8, 9, 10, 10};
    std::vector<size_t> solution_i1 = {3, 4, 5, 6, 7, 8, 9, 10, 11, 11};
    real solution_r = 0.2;
    for (size_t i = 1; i <= 10; i++) {
        auto[i0, i1, r] = SLAdvect::calculate_backward_index(CoordinateAxis::X, i, epsilon, trace_back);
        EXPECT_EQ(i0, solution_i0[i - 1]);
        EXPECT_EQ(i1, solution_i1[i - 1]);
        EXPECT_DOUBLE_EQ(r, solution_r);
    }

    std::vector<size_t> solution_j0 = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11};
    std::vector<size_t> solution_j1 = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12};
    real solution_s = 0.2;
    for (size_t j = 1; j <= 11; j++) {
        auto[j0, j1, s] = SLAdvect::calculate_backward_index(CoordinateAxis::Y, j, epsilon, trace_back);
        EXPECT_EQ(j0, solution_j0[j - 1]);
        EXPECT_EQ(j1, solution_j1[j - 1]);
        EXPECT_DOUBLE_EQ(s, solution_s);
    }

    std::vector<size_t> solution_k0 = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12};
    std::vector<size_t> solution_k1 = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13};
    real solution_t = 0.2;
    for (size_t k = 1; k <= 12; k++) {
        auto[k0, k1, t] = SLAdvect::calculate_backward_index(CoordinateAxis::Z, k, epsilon, trace_back);
        EXPECT_EQ(k0, solution_k0[k - 1]);
        EXPECT_EQ(k1, solution_k1[k - 1]);
        EXPECT_DOUBLE_EQ(t, solution_t);
    }
}

