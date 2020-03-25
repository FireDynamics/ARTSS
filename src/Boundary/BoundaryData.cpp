/// \file 		BoundaryData.cpp
/// \brief 		Data class for boundary data
/// \date 		Oct 08, 2020
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.
#include <iostream>
#include "BoundaryData.h"
#include "../Domain.h"
#include "../Utility/Utility.h"
#include "../Utility/GlobalMacrosTypes.h"

inline static const std::vector<std::string> FieldTypeNames = {"rho", "u", "v", "w", "p", "T"};
inline static const std::vector<std::string> PatchNames = {"front", "back", "bottom", "top", "left", "right"};
inline static const std::vector<std::string> BoundaryConditionNames = {"neumann", "dirichlet", "periodic"};

BoundaryData::BoundaryData() {
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
    for (size_t i = 0; i < numberOfPatches; i++) {
        std::cout << "\t Patch " << getPatchName(static_cast<Patch>(i)) << " with " << getBoundaryConditionName(m_boundaryConditions[i]) << " " << m_values[i] << std::endl;
    }
    //TODO print boundaries
}

//====================================== Matches =================================
// *******************************************************************************
/// \brief  matches string to FieldTypeNames
/// \param	s			string to be matched
// *******************************************************************************
FieldType BoundaryData::matchField(const std::string &s) {

    for (auto fn = 0; fn < FieldTypeNames.size(); fn++) {
        if (FieldTypeNames[fn] == s) return (FieldType) fn;
    }

    std::cout << "unknown field found -> exit" << std::endl;
    std::flush(std::cout);
    std::exit(1);
}

// *******************************************************************************
/// \brief  matches string to PatchNames
/// \param	s			string to be matched
// *******************************************************************************
Patch BoundaryData::matchPatch(const std::string &s) {
    for (auto pn = 0; pn < PatchNames.size(); pn++) {
        if (PatchNames[pn] == s) return (Patch) pn;
    }

    std::cout << "unknown patch found -> exit" << std::endl;
    std::flush(std::cout);
    std::exit(1);
}

// *******************************************************************************
/// \brief  matches string to BoundaryConditionNames
/// \param	s			string to be matched
// *******************************************************************************
BoundaryCondition BoundaryData::matchBoundaryCondition(const std::string &s) {
    for (auto tn = 0; tn < BoundaryConditionNames.size(); tn++) {
        if (BoundaryConditionNames[tn] == s) return (BoundaryCondition) tn;
    }

    std::cout << "unknown BC type found -> exit" << std::endl;
    std::flush(std::cout);
    std::exit(1);
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
/// \param	patches Corresponding patches of boundary condition
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
