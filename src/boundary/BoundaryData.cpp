/// \file       BoundaryData.cpp
/// \brief      Data class for boundary data
/// \date       Oct 08, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.
#include "BoundaryData.h"
#include "../DomainData.h"
#include "../utility/Utility.h"

inline static const std::vector<std::string> boundary_condition_names = {"neumann", "dirichlet", "periodic"};

BoundaryData::BoundaryData() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_values = new real[number_of_patches];
    m_boundary_conditions = new BoundaryCondition[number_of_patches];
}

BoundaryData::~BoundaryData() {
    delete[] m_values;
    delete[] m_boundary_conditions;
}

//====================================== Print =====================================================
// *************************************************************************************************
/// \brief  Print boundary infos
// *************************************************************************************************
void BoundaryData::print() {
#ifndef BENCHMARKING
    for (size_t i = 0; i < number_of_patches; i++) {
        std::string p = PatchObject::get_patch_name(static_cast<Patch>(i));
        std::string bc = get_boundary_condition_name(m_boundary_conditions[i]);
        real val = m_values[i];
        m_logger->info("\t Patch {} with {} {}", p, bc, val);
    }
#endif
}

// *************************************************************************************************
/// \brief  matches string to boundary_condition_names
/// \param  string           string to be matched
// *************************************************************************************************
BoundaryCondition BoundaryData::match_boundary_condition(const std::string &string) {
    for (size_t tn = 0; tn < boundary_condition_names.size(); tn++) {
        if (boundary_condition_names[tn] == string) return static_cast<BoundaryCondition>(tn);
    }
    return UNKNOWN_CONDITION;
}

std::string BoundaryData::get_boundary_condition_name(BoundaryCondition bc) {
    return boundary_condition_names[bc];
}

//=============================== Add boundary condition============================================
// *************************************************************************************************
/// \brief  Add boundary condition
/// \param  patches Corresponding patches of boundary condition
/// \param value boundary condition value
/// \param boundary_condition boundary condition
// *************************************************************************************************
void BoundaryData::add_boundary_condition(
        const std::vector<Patch> &patches,
        real value,
        BoundaryCondition boundary_condition) {
    if (!patches.empty()) {
        m_has_values = true;
    }
    for (Patch patch : patches) {
        size_t p = patch;
        *(m_values + p) = value;
        *(m_boundary_conditions + p) = boundary_condition;
    }
}

