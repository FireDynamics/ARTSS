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
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    auto solver = settings.get("solver/description");
    if (solver.find("NS") != std::string::npos || solver.find("Pressure") != std::string::npos) {
        m_levels = settings.get_size_t("solver/pressure/n_level");
    }
    number_of_inner_cells = new Coordinate<size_t>[m_levels + 1];

    for (size_t axis = 0; axis < number_of_axes; axis++) {
        std::string axis_name = Mapping::get_axis_name(CoordinateAxis(axis));
        start_coords_PD[axis] = settings.get_real("domain_parameters/" + axis_name + "1");
        end_coords_PD[axis] = settings.get_real("domain_parameters/" + axis_name + "2");
        length_PD[axis] = fabs(end_coords_PD[axis] - start_coords_PD[axis]);

        std::string axis_name_lower = axis_name;
        std::transform(axis_name.begin(), axis_name.end(), axis_name_lower.begin(), ::tolower);
        number_of_inner_cells[0][axis] = settings.get_size_t("domain_parameters/n" + axis_name_lower);
    }

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

    calc_MG_values();
#ifndef BENCHMARKING
    print_details();
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
/// \params start_coord_p Start coordinate of physical domain
/// \params end_coord_p End coordinate of physical domain
/// \params old_coord Old coordinate of computational domain
/// \params cell_width  dx,dy,dz
/// \params new_coord New coordinate of computational domain
/// \return bool True if resize is valid
// ***************************************************************************************
bool DomainData::set_new_value(long shift, real start_coord_p, real end_coord_p, real old_coord, real cell_width, real *new_coord) {
    bool update = false;
    if (shift != 0) {
        *new_coord = calc_new_coord(old_coord, shift, cell_width);
        if (*new_coord >= start_coord_p && *new_coord <= end_coord_p) {
            update = true;
        }
    }
    return update;
}

void DomainData::print() {
#ifndef BENCHMARKING
    m_logger->info("-- Domain");
    m_logger->info("Domain size inner cells: {}", *number_of_inner_cells);
    m_logger->info("step size (x|y|z): ({}|{}|{})", get_spacing(CoordinateAxis::X),
                                                    get_spacing(CoordinateAxis::Y),
                                                    get_spacing(CoordinateAxis::Z));
#endif
}

void DomainData::print_details() {
#ifndef BENCHMARKING
    m_logger->debug("############### Domain Data Parameter ###############");
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} (Nx|Ny|Nz): ({}|{}|{})", level,
                get_number_of_cells(CoordinateAxis::X, level),
                get_number_of_cells(CoordinateAxis::Y, level),
                get_number_of_cells(CoordinateAxis::Z, level));
    }
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} (nx|ny|nz): {}", level, number_of_inner_cells[level]);
    }
    m_logger->debug("start coordinates in physical domain (X1|Y1|Z1) {}", start_coords_PD);
    m_logger->debug("end coordinates in physical domain (X2|Y2|Z2) {}", end_coords_PD);

    m_logger->debug("start coordinates in computational domain (x1|y1|z1) {}", start_coords_CD);
    m_logger->debug("end coordinates in computational domain (x1|y1|z1) {}", end_coords_CD);

    m_logger->debug("length of physical domain (Lx|Ly|Lz): {}", length_PD);
    m_logger->debug("length of computational domain (lx|ly|lz): ({}|{}|{})",
                    get_length_CD(CoordinateAxis::X),
                    get_length_CD(CoordinateAxis::Y),
                    get_length_CD(CoordinateAxis::Z));

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} spacing (dx|dy|dz): ({}|{}|{})", level,
                        get_spacing(CoordinateAxis::X, level),
                        get_spacing(CoordinateAxis::Y, level),
                        get_spacing(CoordinateAxis::Z, level));
    }
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} start index computational domain (x|y|z): ({}|{}|{})", level,
                        get_start_index_CD(CoordinateAxis::X, level),
                        get_start_index_CD(CoordinateAxis::Y, level),
                        get_start_index_CD(CoordinateAxis::Z, level));
        m_logger->debug("For Level {} end index computational domain (x|y|z): ({}|{}|{})", level,
                        get_end_index_CD(CoordinateAxis::X, level),
                        get_end_index_CD(CoordinateAxis::Y, level),
                        get_end_index_CD(CoordinateAxis::Z, level));
    }
    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} domain size: {}", level, get_size(level));
    }

    m_logger->debug("--------------- Domain Data Parameter end ---------------");
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

DomainData::~DomainData() {
    delete[] number_of_inner_cells;
}
