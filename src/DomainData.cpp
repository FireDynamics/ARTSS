/// \file       Domain.cpp
/// \brief      XML Domain parameters to variables
/// \date       Jul 16, 2018
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "DomainData.h"

#include <string>


DomainData *DomainData::single = nullptr;  // Singleton

DomainData::DomainData(Settings::Settings const &settings) :
        length_PD(),
        start_coords_PD(),
        end_coords_PD() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(settings, typeid(this).name());
#endif
    auto solver = settings.get("solver/description");
    if (solver.find("NS") != std::string::npos || solver.find("Pressure") != std::string::npos) {
        m_levels = settings.get_size_t("solver/pressure/n_level");
    }
    number_of_inner_cells = new Coordinate<size_t>[m_levels + 1];
    size_t nx = settings.get_size_t("domain_parameters/nx");
    size_t ny = settings.get_size_t("domain_parameters/ny");
    size_t nz = settings.get_size_t("domain_parameters/nz");
    number_of_inner_cells[0].set_coordinate(nx, ny, nz);

    real X1 = settings.get_real("domain_parameters/X1");
    real X2 = settings.get_real("domain_parameters/X2");
    real Y1 = settings.get_real("domain_parameters/Y1");
    real Y2 = settings.get_real("domain_parameters/Y2");
    real Z1 = settings.get_real("domain_parameters/Z1");
    real Z2 = settings.get_real("domain_parameters/Z2");
    start_coords_PD.set_coordinate(X1, Y1, Z1);
    end_coords_PD.set_coordinate(X2, Y2, Z2);

    bool has_computational_domain = settings.get_bool("domain_parameters/enable_computational_domain");
    if (has_computational_domain) {
        real x1 = settings.get_real("domain_parameters/x1");
        real x2 = settings.get_real("domain_parameters/x2");
        real y1 = settings.get_real("domain_parameters/y1");
        real y2 = settings.get_real("domain_parameters/y2");
        real z1 = settings.get_real("domain_parameters/z1");
        real z2 = settings.get_real("domain_parameters/z2");
        start_coords_CD.set_coordinate(x1, y1, z1);
        end_coords_CD.set_coordinate(x2, y2, z2);
    } else {
        start_coords_CD.copy(start_coords_PD);
        end_coords_CD.copy(end_coords_PD);
    }
    for (size_t axis = 0; axis < number_of_axes; axis++) {
        length_PD[axis] = fabs(end_coords_PD[axis] - start_coords_PD[axis]);
    }

    calc_MG_values();
#ifndef BENCHMARKING
    printDetails();
    control();
#endif
}

// =============================== Calculation of MultiGrid Values ========================
// ***************************************************************************************
/// \brief  Calculates amount of Cells in XYZ direction for each multigrid level
// ***************************************************************************************
void DomainData::calc_MG_values() {
    for (size_t level = 1; level < m_levels + 1; ++level) {
        for (size_t axis = 0; axis < number_of_axes; axis++) {
            number_of_inner_cells[level][axis] = (number_of_inner_cells[level - 1][axis] == 1) ? 1 : static_cast<size_t> (std::round(number_of_inner_cells[level - 1][axis] / 2));
        }
    }
}

DomainData *DomainData::getInstance(Settings::Settings const &settings) {
    if (single == nullptr) {
        single = new DomainData(settings);
    }
    return single;
}

size_t DomainData::get_size(size_t level) const {
    size_t size = 1;
    for (size_t axis = 0; axis < number_of_axes; axis++) {
        size *= get_number_of_cells(CoordinateAxis(axis), level);
    }
    return size;
}

// =============================== Resize computational Domain ========================
// ***************************************************************************************
/// \brief  calculates the size of the new calculation domain
/// \params shift_start  Amount of cells which will be resized (left, bottom, front)
/// \params shift_end  Amount of cells which will be resized (right, top, back)
/// \return bool True if computational domain has been resized, False if not
// ***************************************************************************************
bool DomainData::resize(const Coordinate<long> &shift_start, const Coordinate<long> &shift_end) {
    bool update = false;
    real tmp;
    for (size_t axis = 0; axis < number_of_axes; axis++) {
        if (set_new_value(shift_start[axis], start_coords_PD[axis], end_coords_PD[axis], start_coords_CD[axis], get_spacing(CoordinateAxis(axis)), &tmp)) {
            start_coords_CD[axis] = tmp;
            update = true;
        }
        if (set_new_value(shift_end[axis], start_coords_PD[axis], end_coords_PD[axis], end_coords_CD[axis], get_spacing(CoordinateAxis(axis)), &tmp)) {
            end_coords_CD[axis] = tmp;
            update = true;
        }
    }
#pragma acc wait
    if (update) {
#ifndef BENCHMARKING
        m_logger->info("Resize domain at start {} and end {}", shift_start, shift_end);
#endif
        for (size_t axis = 0; axis < number_of_axes; axis++) {
            number_of_inner_cells[0][axis] = static_cast<size_t> (std::round(get_length_CD(CoordinateAxis(axis)) / get_spacing(CoordinateAxis(axis)) + 2));
        }
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
real DomainData::calc_new_coord(real oldCoord, long shift, real cell_width) {
    return oldCoord + static_cast<real>(shift) * cell_width;
}

// =============================== Setting coordinates of new computational domain  ========================
// ***************************************************************************************
/// \brief If resizing is valid set new coordinates
/// \params shift  Amount of cells which will be resized
/// \params startCoord_p Start coordinate of physical domain
/// \params endCoord_p End coordinate of physical domain
/// \params oldCoord Old coordinate of computational domain
/// \params cellWidth  dx,dy,dz
/// \params newCoord New coordinate of computational domain
/// \return bool True if resize is valid
// ***************************************************************************************
bool DomainData::set_new_value(long shift, real startCoord_p, real endCoord_p, real oldCoord, real cell_width, real *newCoord) {
    bool update = false;
    if (shift != 0) {
        *newCoord = calc_new_coord(oldCoord, shift, cell_width);
        if (*newCoord >= startCoord_p && *newCoord <= endCoord_p) {
            update = true;
        }
    }
    return update;
}

void DomainData::print() {
#ifndef BENCHMARKING
    m_logger->info("-- Domain");
    m_logger->info("Domain size inner cells: ({}|{}|{})", get_nx(),
                                                          get_ny(),
                                                          get_nz());
    m_logger->info("step size (x|y|z): ({}|{}|{})", get_dx(),
                                                    get_dy(),
                                                    get_dz());
#endif
}

void DomainData::printDetails() {
#ifndef BENCHMARKING
    m_logger->debug("############### Domain Parameter ###############");
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} Nx: {}, Ny: {}, Nz: {}", level,
                get_Nx(level), get_Ny(level), get_Nz(level));
    }
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} nx: {}, ny {}, nz {}", level,
                get_nx(level), get_ny(level), get_nz(level));
    }

    m_logger->debug("X: ({}|{}) x: ({}|{})", get_X1(), get_X2(),
                                            get_x1(), get_x2());
    m_logger->debug("Y: ({}|{}) y: ({}|{})", get_Y1(), get_Y2(),
                                            get_y1(), get_y2());
    m_logger->debug("Z: ({}|{}) z: ({}|{})", get_Z1(), get_Z2(),
                                            get_z1(), get_z2());

    m_logger->debug("Lx: {}, Ly: {}, Lz: {}", get_Lx(), get_Ly(), get_Lz());
    m_logger->debug("lx: {}, ly: {}, lz: {}", get_lx(), get_ly(), get_lz());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} dx: {}, dy: {}, dz: {}", level,
                get_dx(level), get_dy(level), get_dz(level));
    }
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("for Level {} X: ({}|{}) Y: ({}|{}) Z: ({}|{})", level,
                get_index_x1(level), get_index_x2(level),
                get_index_y1(level), get_index_y2(level),
                get_index_z1(level), get_index_z2(level));
    }
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug(" For Level {} domain size: {}", level, get_size(level));
    }

    m_logger->debug("--------------- Domain Parameter end ---------------");
#endif
}

void DomainData::control() {
#ifndef BENCHMARKING
    if (m_levels > 0) {
        int minimum_amount_of_cells = static_cast<int>(std::pow(2, m_levels));
        for (size_t axis = 0; axis < number_of_axes; axis++) {
            size_t number_of_cells = get_number_of_inner_cells(CoordinateAxis(axis));
            if (number_of_cells % minimum_amount_of_cells != 0) {
                std::string spacing_name = axis_names[axis];
                std::transform(spacing_name.begin(), spacing_name.end(), spacing_name.begin(), ::tolower);
                m_logger->warn("n{} ({}) has to be a multiple of 2^levels = {}. "
                               "Consider changing n{} or levels of multigrid solver. "
                               "Otherwise unexpected behaviour may occur.",
                               spacing_name, number_of_cells, minimum_amount_of_cells,
                               spacing_name);
            }
        }
    }
#endif
}
