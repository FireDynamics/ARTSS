/// \file       BoundaryData.cpp
/// \brief      Data class for boundary data
/// \date       Oct 08, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.
#include "BoundaryData.h"
#include "../Domain.h"
#include "../utility/Utility.h"

inline static const std::vector<std::string> FieldTypeNames = {"rho", "u", "v", "w", "p", "T"};
inline static const std::vector<std::string> PatchNames = {"front", "back", "bottom", "top", "left", "right"};
inline static const std::vector<std::string> BoundaryConditionNames = {"neumann", "dirichlet", "periodic"};

BoundaryData::BoundaryData() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_values = new real[numberOfPatches];
    m_boundaryConditions = new BoundaryCondition[numberOfPatches];
}

BoundaryData::~BoundaryData() {
    delete[] m_values;
    delete[] m_boundaryConditions;
}

//====================================== Print =================================
// *******************************************************************************
/// \brief  Print boundary infos
// *******************************************************************************
void BoundaryData::print() {
#ifdef BENCHMARKING
    return;
#else
    for (size_t i = 0; i < numberOfPatches; i++) {
        std::string p = getPatchName(static_cast<Patch>(i));
        std::string bc = getBoundaryConditionName(m_boundaryConditions[i]);
        real val = m_values[i];
        m_logger->info("\t Patch {} with {} {}", p , bc, val);
    }
#endif
}

//====================================== Matches =================================
// *******************************************************************************
/// \brief  matches string to FieldTypeNames
/// \param  s           string to be matched
// *******************************************************************************
FieldType BoundaryData::matchField(const std::string &s) {
    for (size_t fn = 0; fn < FieldTypeNames.size(); fn++) {
        if (FieldTypeNames[fn] == s) return (FieldType) fn;
    }

    // unknown fieldtype => die
    return UNKOWN_FIELD;
}

// *******************************************************************************
/// \brief  matches string to PatchNames
/// \param  s           string to be matched
// *******************************************************************************
Patch BoundaryData::matchPatch(const std::string &s) {
    for (size_t pn = 0; pn < PatchNames.size(); pn++) {
        if (PatchNames[pn] == s) return (Patch) pn;
    }

    return UNKOWN_PATCH;
}

// *******************************************************************************
/// \brief  matches string to BoundaryConditionNames
/// \param  s           string to be matched
// *******************************************************************************
BoundaryCondition BoundaryData::matchBoundaryCondition(const std::string &s) {
    for (size_t tn = 0; tn < BoundaryConditionNames.size(); tn++) {
        if (BoundaryConditionNames[tn] == s) return (BoundaryCondition) tn;
    }

    return UNKOWN_CONDITION;
}

std::string BoundaryData::getFieldTypeName(FieldType f) {
    return FieldTypeNames[f];
}

std::string BoundaryData::getBoundaryConditionName(BoundaryCondition bc) {
    return BoundaryConditionNames[bc];
}

std::string BoundaryData::getPatchName(Patch p) {
    return PatchNames[p];
}

//=============================== Add boundary condition=========================
// *******************************************************************************
/// \brief  Add boundary condition
/// \param  patches Corresponding patches of boundary condition
/// \param value boundary condition value
/// \param boudnaryCondition boundary condition
// *******************************************************************************
void BoundaryData::addBoundaryCondition(const std::vector<Patch> &patches, real value, BoundaryCondition boundaryCondition) {
    if (!patches.empty()) {
        m_hasValues = true;
    }
    for (Patch patch: patches) {
        size_t p = patch;
        *(m_values + p) = value;
        *(m_boundaryConditions + p) = boundaryCondition;
    }
}
