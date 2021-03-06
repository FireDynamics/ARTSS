/// \file       BoundaryDataController.h
/// \brief      Controller class for boundary data
/// \date       Dec 09, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "BoundaryDataController.h"
#include <string>
#include "../utility/Utility.h"
#include "../boundaryCondition/DomainBoundary.h"
#include "../boundaryCondition/ObstacleBoundary.h"

BoundaryDataController::BoundaryDataController() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_boundary_data = new BoundaryData *[numberOfFieldTypes];
    for (size_t i = 0; i < numberOfFieldTypes; i++) {
        *(m_boundary_data + i) = new BoundaryData();
    }
}

BoundaryDataController::~BoundaryDataController() {
    for (size_t i = 0; i < numberOfFieldTypes; i++) {
        delete m_boundary_data[i];
    }
    delete[] m_boundary_data;
}

// =================================== Add boundary data ===========================================
// *************************************************************************************************
/// \brief  Parses boundary data of XML tree to boundary data object
/// \param  xml_element Pointer to XML element
// *************************************************************************************************
void BoundaryDataController::add_boundary_data(tinyxml2::XMLElement *xml_element) {
    std::vector<std::string> fieldStrings = Utility::split(xml_element->Attribute("field"), ',');
    std::vector<std::string> patchStrings = Utility::split(xml_element->Attribute("patch"), ',');
    std::vector<Patch> patches;
    patches.reserve(patchStrings.size());

    for (const std::string &p : patchStrings) {
        patches.push_back(BoundaryData::match_patch(p));
    }

    BoundaryCondition boundaryCondition = BoundaryData::match_boundary_condition(xml_element->Attribute("type"));
    auto value = xml_element->DoubleAttribute("value");

    for (const std::string &f : fieldStrings) {
        FieldType fieldType = BoundaryData::match_field(f);
        m_boundary_data[fieldType]->add_boundary_condition(patches, value, boundaryCondition);
    }
}

// ============================================== Print ============================================
// *************************************************************************************************
/// \brief  Prints info of boundary data
// *************************************************************************************************
void BoundaryDataController::print() {
#ifndef BENCHMARKING
    for (size_t i = 0; i < numberOfFieldTypes; i++) {
        auto boundary = *(m_boundary_data + i);
        if (!boundary->is_empty()) {
            m_logger->info("--- found boundary conditions for field {} ({}): ",
                           BoundaryData::get_field_type_name(static_cast<FieldType>(i)), i);
            boundary->print();
        }
    }
#endif
}

//======================================== Apply boundary condition ================================
// *************************************************************************************************
/// \brief  Applies boundary condition for domain boundary if the field is needed
/// \param  dataField     Field
/// \param  index_fields  List of indices for each patch
/// \param  patch_starts  List of start indices
/// \param  patch_ends    List of end indices
/// \param  field_type    Type of field
/// \param  level         Multigrid level
/// \param  sync          synchronous kernel launching (true, default: false)
// *************************************************************************************************
void BoundaryDataController::apply_boundary_condition(
        real *data,
        size_t **index_fields,
        size_t *patch_start, size_t *patch_end,
        FieldType field_type,
        size_t level,
        bool sync) {
    if (!((static_cast<BoundaryData *> (*(m_boundary_data + field_type)))->is_empty())) {
        DomainBoundary::apply_boundary_condition(
                data,
                index_fields,
                patch_start, patch_end,
                level,
                m_boundary_data[field_type],
                sync);
    }
}

//=========================== Apply obstacle boundary condition ====================================
// *************************************************************************************************
/// \brief  Applies boundary condition for obstacle boundary if the field is needed
/// \param  dataField     Field
/// \param  index_fields  List of indices for each patch
/// \param  patch_starts  List of start indices
/// \param  patch_ends    List of end indices
/// \param  field_type    Type of field
/// \param  level         Multigrid level
/// \param  id            ID of obstacle
/// \param  sync          synchronous kernel launching (true, default: false)
// *************************************************************************************************
void BoundaryDataController::apply_boundary_condition_obstacle(
        real *data,
        size_t **index_fields,
        size_t *patch_start, size_t *patch_end,
        FieldType field_type,
        size_t level,
        size_t id,
        bool sync) {
    if (!(static_cast<BoundaryData *> (*(m_boundary_data + field_type)))->is_empty()) {
#ifndef BENCHMARKING
        m_logger->debug("apply obstacle boundary conditions of {}",
                        BoundaryData::get_field_type_name(static_cast<FieldType>(field_type)));
#endif
        ObstacleBoundary::apply_boundary_condition(
                data,
                index_fields,
                patch_start, patch_end,
                level,
                m_boundary_data[field_type],
                id,
                sync);
    }
}

//==================================== get used fields =============================================
// *************************************************************************************************
/// \brief returns only the fields which are use. For a pure advection case this function will only
/// return the velocity field
/// \return returns a vector with the used fields
// *************************************************************************************************
std::vector<FieldType> BoundaryDataController::get_used_fields() {
    std::vector<FieldType> v_fields;
    v_fields.reserve(numberOfFieldTypes);
    for (size_t fieldType = 0; fieldType < numberOfFieldTypes; fieldType++) {
        if (!(static_cast<BoundaryData *> (*(m_boundary_data + fieldType)))->is_empty()) {
            v_fields.push_back(static_cast<FieldType>(fieldType));
        }
    }
    return v_fields;
}
