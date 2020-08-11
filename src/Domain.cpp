/// \file       Domain.cpp
/// \brief      XML Domain parameters to variables
/// \date       Jul 16, 2018
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Domain.h"

Domain *Domain::single = nullptr; //Singleton

Domain::Domain() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    auto params = Parameters::getInstance();
    auto solver = params->get("solver/description");
    if (solver.find("NS") != std::string::npos || solver.find("Pressure") != std::string::npos) {
        m_levels = static_cast<size_t> (params->get_int("solver/pressure/n_level"));
    }
    m_nx = new size_t[m_levels + 1];
    m_ny = new size_t[m_levels + 1];
    m_nz = new size_t[m_levels + 1];

    m_nx[0] = static_cast<size_t> (params->get_int("domain_parameters/nx")) + 2;
    m_ny[0] = static_cast<size_t> (params->get_int("domain_parameters/ny")) + 2;
    m_nz[0] = static_cast<size_t> (params->get_int("domain_parameters/nz")) + 2;

    m_x1 = params->get_real("domain_parameters/x1");
    m_x2 = params->get_real("domain_parameters/x2");
    m_y1 = params->get_real("domain_parameters/y1");
    m_y2 = params->get_real("domain_parameters/y2");
    m_z1 = params->get_real("domain_parameters/z1");
    m_z2 = params->get_real("domain_parameters/z2");
    m_X1 = params->get_real("domain_parameters/X1");
    m_X2 = params->get_real("domain_parameters/X2");
    m_Y1 = params->get_real("domain_parameters/Y1");
    m_Y2 = params->get_real("domain_parameters/Y2");
    m_Z1 = params->get_real("domain_parameters/Z1");
    m_Z2 = params->get_real("domain_parameters/Z2");

    calc_MG_values();
#ifndef BENCHMARKING
    //printDetails();
#endif
}

// =============================== Calculation of MultiGrid Values ========================
// ***************************************************************************************
/// \brief  Calculates amount of Cells in XYZ direction for each multigrid level
// ***************************************************************************************
void Domain::calc_MG_values() {

    for (size_t l = 1; l < m_levels + 1; ++l) {
        m_nx[l] = (m_nx[l - 1] == 3) ? 3 : static_cast<size_t> (std::round((m_nx[l - 1] - 2) / 2 + 2));
        m_ny[l] = (m_ny[l - 1] == 3) ? 3 : static_cast<size_t> (std::round((m_ny[l - 1] - 2) / 2 + 2));
        m_nz[l] = (m_nz[l - 1] == 3) ? 3 : static_cast<size_t> (std::round((m_nz[l - 1] - 2) / 2 + 2));
    }
}

Domain *Domain::getInstance() {
    if (single == nullptr) {
        single = new Domain();
    }
    return single;
}

// =============================== Resize computational Domain ========================
// ***************************************************************************************
/// \brief  calculates the size of the new calculation domain
/// \params shift_x1  Amount of cells which will be resized in x direction
/// \params shift_x2  Amount of cells which will be resized in x direction
/// \params shift_y1  Amount of cells which will be resized in y direction
/// \params shift_y2  Amount of cells which will be resized in y direction
/// \params shift_z1  Amount of cells which will be resized in z direction
/// \params shift_z2  Amount of cells which will be resized in z direction
/// \return bool True if computational domain has been resized, False if not
// ***************************************************************************************
bool Domain::resize(long shift_x1, long shift_x2, long shift_y1, long shift_y2, long shift_z1, long shift_z2) {
    auto dx = get_dx();
    auto dy = get_dy();
    auto dz = get_dz();

    bool update = false;
    real tmp;
    if (set_new_value(shift_x1, m_X1, m_X2, m_x1, dx, &tmp)) {
        m_x1 = tmp;
        update = true;
    }
    if (set_new_value(shift_x2, m_X1, m_X2, m_x2, dx, &tmp)) {
        m_x2 = tmp;
        update = true;
    }
    if (set_new_value(shift_y1, m_Y1, m_Y2, m_y1, dy, &tmp)) {
        m_y1 = tmp;
        update = true;
    }
    if (set_new_value(shift_y2, m_Y1, m_Y2, m_y2, dy, &tmp)) {
        m_y2 = tmp;
        update = true;
    }
    if (set_new_value(shift_z1, m_Z1, m_Z2, m_z1, dz, &tmp)) {
        m_z1 = tmp;
        update = true;
    }
    if (set_new_value(shift_z2, m_Z1, m_Z2, m_z2, dz, &tmp)) {
        m_z2 = tmp;
        update = true;
    }
#pragma acc wait
    if (update) {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(Domain).name());
        m_logger->info("Resize domain: {}|{} {}|{} {}|{}", shift_x1, shift_x2,
                                                         shift_y1, shift_y2,
                                                         shift_z1, shift_z2);
#endif
        m_nx[0] = static_cast<size_t> (std::round(get_lx() / dx + 2));
        m_ny[0] = static_cast<size_t> (std::round(get_ly() / dy + 2));
        m_nz[0] = static_cast<size_t> (std::round(get_lz() / dz + 2));

        calc_MG_values();
    }
    return update;
}

// =============================== Calculation of new Coordinates  ========================
// ***************************************************************************************
/// \brief Calculates new coordinates depending on shift variable
/// \params oldCoord  coordinate of old computational domain
/// \params shift  Amount of cells which will be resized
/// \params cellWidth  dx,dy,dz
/// \return coordinate of new computational domain
// ***************************************************************************************
real Domain::calc_new_coord(real oldCoord, long shift, real cell_width) {
    return oldCoord + static_cast<real>(shift) * cell_width;
}

// =============================== Setting coordinates of new computational domain  ========================
// ***************************************************************************************
/// \brief If resizing is valid set new coordinates
/// \params shift  Amount of cells which will be resized
/// \params startCoord_p Start coordinate of physical domain
/// \params endCoord_p End coordinate of physical domain
/// \params oldCoord Old coordiante of computational domain
/// \params cellWidth  dx,dy,dz
/// \params newCoord New coordiante of computational domain
/// \return bool True if resize is valid
// ***************************************************************************************
bool Domain::set_new_value(long shift, real startCoord_p, real endCoord_p, real oldCoord, real cell_width, real *newCoord) {
    bool update = false;
    if (shift != 0) {
        *newCoord = calc_new_coord(oldCoord, shift, cell_width);
        if (*newCoord >= startCoord_p && *newCoord <= endCoord_p) {
            update = true;
        }
    }
    return update;
}

void Domain::print() {
#ifdef BENCHMARKING
    return;
#else
    m_logger->info("-- Domain");
    m_logger->info("Domain size inner cells: ({}|{}{})", get_nx() - 2,
                                                         get_ny() - 2,
                                                         get_nz() - 2);
    m_logger->info("step size (x|y|z): ({}|{}|{})", get_dx(),
                                                    get_dy(),
                                                    get_dz());
#endif
}

void Domain::printDetails() {
#ifdef BENCHMARKING
    return;
#else
    m_logger->info("############### Domain Parameter ###############");
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->info("For Level {} Nx: {}, Ny: {}, Nz: {}", level,
                get_Nx(level), get_Ny(level), get_Nz(level));
    }
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->info("For Level {} nx: {}, ny {}, nz {}", level,
                get_nx(level), get_ny(level), get_nz(level));
    }

    m_logger->info("X: ({}|{}) x: ({}|{})", get_X1(), get_X2(),
                                            get_x1(), get_x2());
    m_logger->info("Y: ({}|{}) y: ({}|{})", get_Y1(), get_Y2(),
                                            get_y1(), get_y2());
    m_logger->info("Z: ({}|{}) z: ({}|{})", get_Z1(), get_Z2(),
                                            get_z1(), get_z2());

    m_logger->info("Lx: {}, Ly: {}, Lz: {}", get_Lx(), get_Ly(), get_Lz());
    m_logger->info("lx: {}, ly: {}, lz: {}", get_lx(), get_ly(), get_lz());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->info("For Level {} dx: {}, dy: {}, dz: {}", level,
                get_dx(level), get_dy(level), get_dz(level));
    }
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->info("for Level {} X: ({}|{}) Y: ({}|{}) Z: ({}|{})", level,
                get_index_x1(level), get_index_x2(level),
                get_index_y1(level), get_index_y2(level),
                get_index_z1(level), get_index_z2(level));
    }
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->info(" For Level {} domain size: {}", level, get_size(level));
    }

    m_logger->info("--------------- Domain Parameter end ---------------");
#endif
}
