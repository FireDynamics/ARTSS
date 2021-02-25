/// \file       GaussFunction.cpp
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "GaussFunction.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"


real constexpr det3(real a1, real a2, real a3,
        real b1, real b2, real b3,
        real c1, real c2, real c3) {
    return a1 * (b2 * c3 - b3 * c2)
        - b1 * (a2 * c3 - a3 * c2)
        + c1 * (a2 * b3 - a3 * b2);
}

real constexpr eps = 10E-6;

constexpr int cuboid_points[8][3] = {
    {0, 2, 4},
    {1, 2, 4},
    {1, 3, 4},
    {0, 3, 4},
    {0, 2, 5},
    {1, 2, 5},
    {1, 3, 5},
    {0, 3, 5}
};

constexpr int surface_points[6][3] = {
    {0, 1, 4},
    {1, 2, 5},
    {3, 2, 6},
    {1, 3, 7},
    {0, 1, 3},
    {4, 5, 7},
};

GaussFunction::GaussFunction(real HRR, real cp, real x0, real y0, real z0, real sigma_x, real sigma_y, real sigma_z, real tau) {
    m_tau = tau;
    m_field_spatial_values = new Field(FieldType::RHO, 0.);
    create_spatial_values(HRR, cp, x0, y0, z0, sigma_x, sigma_y, sigma_z);
}

GaussFunction::~GaussFunction() {
    auto data_spatial = m_field_spatial_values->data;
    size_t size = Domain::getInstance()->get_size();
#pragma acc exit data delete(data_spatial[:size])
    delete m_field_spatial_values;
}

void GaussFunction::update_source(Field *out, real t_cur) {
    size_t size = Domain::getInstance()->get_size();
    auto data_out = out->data;
    auto data_spatial = m_field_spatial_values->data;

#pragma acc data present(data_out[:size], data_spatial[:size])
    {
#pragma acc parallel loop independent present(data_out[:size], data_spatial[:size]) async
        for (size_t i = 0; i < size; i++) {
            data_out[i] = data_spatial[i] * get_time_value(t_cur);
        }
#pragma acc wait
    }
}

// ***************************************************************************************
/// \brief  Volumetric Gaussian temperature source in energy equation
/// \param  out   energy source
/// \param  HRR   total heat release rate
/// \param  cp    heat capacity
/// \param  x0    center of Gaussian (x-direction)
/// \param  y0    center of Gaussian (y-direction)
/// \param  z0    center of Gaussian (z-direction)
/// \param  sigma Radius of Gaussian
// ***************************************************************************************
void GaussFunction::create_spatial_values(real HRR, real cp,
        real x0, real y0, real z0,
        real sigma_x, real sigma_y, real sigma_z) {
    auto domain = Domain::getInstance();
    auto bsize = domain->get_size();
    // local variables and parameters for GPU
    auto d_out = m_field_spatial_values->data;
    auto level = m_field_spatial_values->get_level();

    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);

    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    real dx = domain->get_dx(level);
    real dy = domain->get_dy(level);
    real dz = domain->get_dz(level);

    // get parameters for Gaussian
    real sigma_x_2 = 2 * sigma_x * sigma_x;
    real r_sigma_x_2 = 1. / sigma_x_2;
    real sigma_y_2 = 2 * sigma_y * sigma_y;
    real r_sigma_y_2 = 1. / sigma_y_2;
    real sigma_z_2 = 2 * sigma_z * sigma_z;
    real r_sigma_z_2 = 1. / sigma_z_2;

    // set Gaussian to cells
    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();

    const auto multigrid = boundary->getMultigrid();
    const auto obst_size = multigrid->getSize_NumberOfObstacles();
    const auto obst_list = multigrid->getObstacles();
    auto bsize_i = boundary->getSize_innerList();

    auto i0 = (x0 - X1) / dx;
    auto j0 = (y0 - Y1) / dy;
    auto k0 = (z0 - Z1) / dz;

    real V = 0.;
    for (size_t l = 0; l < bsize_i; ++l) {
        const size_t idx = d_iList[l];
        int k = getCoordinateK(idx, Nx, Ny);
        int j = getCoordinateJ(idx, Nx, Ny, k);
        int i = getCoordinateI(idx, Nx, Ny, j, k);

        real di = (i0 - i);
        real dj = (j0 - j);
        real dk = (k0 - k);

        bool blocked = false;
        for (auto obst_id=0; obst_id < obst_size; ++obst_id) {
            auto obst = obst_list[level][obst_id];
            auto point_i1 = obst->getCoordinates_i1();
            auto point_i2 = obst->getCoordinates_i2();
            auto point_j1 = obst->getCoordinates_j1();
            auto point_j2 = obst->getCoordinates_j2();
            auto point_k1 = obst->getCoordinates_k1();
            auto point_k2 = obst->getCoordinates_k2();
            real indeces[6]  = {static_cast<real>(point_i1), static_cast<real>(point_i2),
                                static_cast<real>(point_j1), static_cast<real>(point_j2),
                                static_cast<real>(point_k1), static_cast<real>(point_k2)};

            for (int surface_id=0; surface_id < 6; ++surface_id) {
                auto s = surface_points[surface_id];
                auto p1 = cuboid_points[s[0]];
                auto p2 = cuboid_points[s[1]];
                auto p3 = cuboid_points[s[2]];

                // surface vector 1
                auto svi1 = indeces[p2[0]] - indeces[p1[0]];  // saving dx
                auto svj1 = indeces[p2[1]] - indeces[p1[1]];
                auto svk1 = indeces[p2[2]] - indeces[p1[2]];

                // surface vector 2
                auto svi2 = indeces[p3[0]] - indeces[p1[0]];
                auto svj2 = indeces[p3[1]] - indeces[p1[1]];
                auto svk2 = indeces[p3[2]] - indeces[p1[2]];

                // p1 + l*sv1 + m*sv2 = n*l + c <=>
                // l*sv1 + m*sv2 - n*l = c - p1
                auto det_A = det3(svi1, svj1, svk1,
                    svi2, svj2, svk2,
                    -di, -dj, -dk);

                // okay det is small, its never gonna meet
                if (fabs(det_A) < eps) {
                    continue;
                }

                // rhs of les (c - p1)
                auto ddi1 = static_cast<real>(i) - indeces[p1[0]];
                auto ddj1 = static_cast<real>(j) - indeces[p1[1]];
                auto ddk1 = static_cast<real>(k) - indeces[p1[2]];

                auto det_Ax = det3(ddi1, ddj1, ddk1,
                    svi2, svj2, svk2,
                    -di, -dj, -dk);
                auto det_Ay = det3(svi1, svj1, svk1,
                    ddi1, ddj1, ddk1,
                    -di, -dj, -dk);
                auto det_Az = det3(svi1, svj1, svk1,
                    svi2, svj2, svk2,
                    ddi1, ddj1, ddk1);

                auto sx = - det_Ax / det_A;  // sx : l
                auto sy = - det_Ay / det_A;  // sy : m
                auto sz = - det_Az / det_A;  // sz : n

                blocked = sx > -eps && sx < 1.0
                        && sy > -eps && sy < 1.0
                        && sz > -eps && sz < 1.0;
                
                // std::cout << svi1 << "," << svi2 << "," << indeces[p1[0]] << std::endl;
                // std::cout << svj1 << "," << svj2 << "," << indeces[p1[1]] << std::endl;
                // std::cout << svk1 << "," << svk2 << "," << indeces[p1[2]] << std::endl;

                // std::cout << surface_id << std::endl;

                // std::cout << sx << std::endl;
                // std::cout << sy << std::endl;
                // std::cout << sz << std::endl;

                if (blocked) {
                    break;
                }
            }
        }

        real x_i = dx * di;
        real y_j = dy * dj;
        real z_k = dz * dk;

        real expr = std::exp(-(r_sigma_x_2 * (x_i * x_i)
                                + r_sigma_y_2 * (y_j * y_j)
                                + r_sigma_z_2 * (z_k * z_k)));

        if (blocked) {
            // std::cout << "XXX" << i << "," << j << "," << k << ",0" << std::endl;
            d_out[idx] = 0.0;
        } else {
            // std::cout << "XXX" << i << "," << j << "," << k << ",1" << std::endl;
            d_out[idx] = expr;
        }

        V += expr * dx * dy * dz;
    }

    const real HRRrV = HRR / V;  // in case of concentration Ys*HRR
    const real factor = HRRrV / cp;
    for (size_t l = 0; l < bsize_i; ++l) {
        const size_t idx = d_iList[l];
        d_out[idx] = factor * d_out[idx];
    }

#pragma acc enter data copyin(d_out[:bsize])
}

// ============================= Ramp up function for HRR source =========================
// ***************************************************************************************
/// \brief  Ramp up function (in time) for Gaussian source in energy equation
/// \param  t time
// ***************************************************************************************
real GaussFunction::get_time_value(real t_cur) {
    return tanh(t_cur / m_tau);
}
