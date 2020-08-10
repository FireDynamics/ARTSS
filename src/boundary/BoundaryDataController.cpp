/// \file 		BoundaryDataController.h
/// \brief 		Controll class for boundary data
/// \date 		Dec 09, 2019
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "BoundaryDataController.h"
#include "../utility/Utility.h"
#include "../boundaryCondition/DomainBoundary.h"
#include "../boundaryCondition/ObstacleBoundary.h"

BoundaryDataController::BoundaryDataController() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_boundaryData = new BoundaryData *[numberOfFieldTypes];
    for (size_t i = 0; i < numberOfFieldTypes; i++) {
        *(m_boundaryData + i) = new BoundaryData();
    }
}

BoundaryDataController::~BoundaryDataController() {
    for (size_t i = 0; i < numberOfFieldTypes; i++) {
        delete m_boundaryData[i];
    }
    delete[] m_boundaryData;
}

// ================================= Add boundary data ===========================================
// ***************************************************************************************
/// \brief  Parses boundary data of XML tree to bounday data object
/// \param 	xmlElement Pointer to XML element
// ***************************************************************************************
void BoundaryDataController::addBoundaryData(tinyxml2::XMLElement *xmlElement) {
    std::vector<std::string> fieldStrings = Utility::split(xmlElement->Attribute("field"), ',');
    std::vector<std::string> patchStrings = Utility::split(xmlElement->Attribute("patch"), ',');
    std::vector<Patch> patches;
    patches.reserve(patchStrings.size());

    for (const std::string &p: patchStrings) {
        patches.push_back(BoundaryData::matchPatch(p));
    }

    BoundaryCondition boundaryCondition = BoundaryData::matchBoundaryCondition(xmlElement->Attribute("type"));
    auto value = xmlElement->DoubleAttribute("value");

    for (const std::string &f: fieldStrings) {
        FieldType fieldType = BoundaryData::matchField(f);
        m_boundaryData[fieldType]->addBoundaryCondition(patches, value, boundaryCondition);
    }
}

// ================================= Print ===========================================
// ***************************************************************************************
/// \brief  Prints info of boundary data
// ***************************************************************************************
void BoundaryDataController::print() {
#ifdef BENCHMARKING
    return;
#else
    for (size_t i = 0; i < numberOfFieldTypes; i++) {
        auto boundary = *(m_boundaryData + i);
        if (!boundary->isEmpty()) {
            m_logger->info("--- found boundary conditions for field {} ({}): ", BoundaryData::getFieldTypeName(static_cast<FieldType>(i)), i);
            boundary->print();
        }
    }
#endif
}

//======================================== Apply boundary condition ====================================
// ***************************************************************************************
/// \brief  Applies boundary condition for domain boundary if the field is needed
/// \param  dataField	Field
/// \param  indexFields List of indices for each patch
/// \param  patch_starts List of start indices
/// \param  patch_ends List of end indices
/// \param  fieldType Type of field
/// \param  level Multigrid level
/// \param  sync synchronous kernel launching (true, default: false)
// ***************************************************************************************
void BoundaryDataController::applyBoundaryCondition(real *data, size_t **indexFields, size_t *patch_start, size_t *patch_end, FieldType fieldType, size_t level, bool sync) {
    if (!((BoundaryData *) *(m_boundaryData + fieldType))->isEmpty()) {
        DomainBoundary::applyBoundaryCondition(data, indexFields, patch_start, patch_end, level, m_boundaryData[fieldType], sync);
    }
}

//=========================== Apply obstacle boundary condition ==========================
// ***************************************************************************************
/// \brief  Applies boundary condition for obstacle boundary if the field is needed
/// \param  dataField	Field
/// \param  indexFields List of indices for each patch
/// \param  patch_starts List of start indices
/// \param  patch_ends List of end indices
/// \param  fieldType Type of field
/// \param  level Multigrid level
/// \param  id ID of obstacle
/// \param  sync synchronous kernel launching (true, default: false)
// ***************************************************************************************
void BoundaryDataController::applyBoundaryConditionObstacle(real *data, size_t **indexFields, size_t *patch_start, size_t *patch_end, FieldType fieldType, size_t level, size_t id, bool sync) {
    if (!((BoundaryData *) *(m_boundaryData + fieldType))->isEmpty()) {
        ObstacleBoundary::applyBoundaryCondition(data, indexFields, patch_start, patch_end, level, m_boundaryData[fieldType], id, sync);
    }
}

std::vector<FieldType> BoundaryDataController::get_used_fields() {
    std::vector<FieldType> v_fields;
    for (size_t fieldType = 0; fieldType < numberOfFieldTypes; fieldType++) {
        if (!((BoundaryData *) *(m_boundaryData + fieldType))->isEmpty()){
            v_fields.push_back(static_cast<FieldType>(fieldType));
        }
    }
    return v_fields;
}
