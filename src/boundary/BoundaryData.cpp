/// \file       BoundaryData.cpp
/// \brief      Data class for boundary data
/// \date       Oct 08, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.
#include "BoundaryData.h"
#include "../Domain.h"
#include "../utility/Utility.h"

inline static const std::vector<std::string> field_type_names = {"rho", "u", "v", "w", "p", "T"};
inline static const std::vector<std::string> patch_names = {"front", "back", "bottom", "top", "left", "right"};
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
        std::string p = get_patch_name(static_cast<Patch>(i));
        std::string bc = get_boundary_condition_name(m_boundary_conditions[i]);
        real val = m_values[i];
        m_logger->info("\t Patch {} with {} {}", p, bc, val);
    }
#endif
}

//====================================== Matches ===================================================
// *************************************************************************************************
/// \brief  matches string to field_type_names
/// \param  string           string to be matched
// *************************************************************************************************
FieldType BoundaryData::match_field(const std::string &string) {
    for (size_t fn = 0; fn < field_type_names.size(); fn++) {
        if (field_type_names[fn] == string) return (FieldType) fn;
    }
    return UNKNOWN_FIELD;
}

// *************************************************************************************************
/// \brief  matches string to patch_names
/// \param  string           string to be matched
// *************************************************************************************************
Patch BoundaryData::match_patch(const std::string &string) {
    for (size_t pn = 0; pn < patch_names.size(); pn++) {
        if (patch_names[pn] == string) return (Patch) pn;
    }
    return UNKNOWN_PATCH;
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

std::string BoundaryData::get_field_type_name(FieldType f) {
    return field_type_names[f];
}

std::string BoundaryData::get_boundary_condition_name(BoundaryCondition bc) {
    return boundary_condition_names[bc];
}

std::string BoundaryData::get_patch_name(Patch p) {
    return patch_names[p];
}

//=============================== Add boundary condition============================================
// *************************************************************************************************
/// \brief  Add boundary condition
/// \param  patches Corresponding patches of boundary condition
/// \param value boundary condition value
/// \param boudnaryCondition boundary condition
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

