/// \file       BoundaryData.cpp
/// \brief      Data class for boundary data
/// \date       Oct 08, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.
#include "BoundaryData.h"
#include "DomainData.h"


BoundaryData::BoundaryData() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_values.resize(number_of_patches);
    m_boundary_conditions.resize(number_of_patches);
}

//====================================== Print =====================================================
// *************************************************************************************************
/// \brief  Print boundary infos
// *************************************************************************************************
void BoundaryData::print() const {
#ifndef BENCHMARKING
    for (size_t i = 0; i < number_of_patches; i++) {
        std::string p = Mapping::get_patch_name(Patch(i));
        std::string bc = Mapping::get_boundary_condition_name(m_boundary_conditions[i]);
        real val = m_values[i];
        m_logger->info("\t Patch {} with {} {}", p, bc, val);
    }
#endif
}

//=============================== Add boundary condition============================================
// *************************************************************************************************
/// \brief  Add boundary condition
/// \param  patches Corresponding patches of boundary condition
/// \param value boundary condition value
/// \param boundary_condition boundary condition
// *************************************************************************************************
void BoundaryData::add_boundary_condition(
        Patch const &patch,
        real value,
        BoundaryCondition const &boundary_condition) {
    m_has_values = true;
    m_values[patch] = value;
    m_boundary_conditions[patch] = boundary_condition;
}

