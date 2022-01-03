/// \file       BoundaryDataController.h
/// \brief      Controller class for boundary data
/// \date       Dec 09, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "BoundaryDataController.h"

#include <string>

#include "../boundaryCondition/DomainBoundary.h"
#include "../boundaryCondition/ObstacleBoundary.h"
#include "../boundaryCondition/SurfaceBoundary.h"


BoundaryDataController::BoundaryDataController(const std::vector<Settings::boundary> &boundary) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_boundary_data.resize(number_of_field_types);
    for (size_t i = 0; i < number_of_field_types; i++) {
        m_boundary_data.emplace_back(BoundaryData());
    }

    for (const auto &b : boundary) {
        add_boundary_data(b);
    }
}

// =================================== Add boundary data ===========================================
// *************************************************************************************************
/// \brief  Parses boundary data of XML tree to boundary data object
/// \param  xml_element Pointer to XML element
// *************************************************************************************************
void BoundaryDataController::add_boundary_data(const Settings::boundary &boundary) {
    for (const auto &field_type : boundary.field_type) {
        for (const auto &patch : boundary.patch) {
            m_boundary_data[field_type].add_boundary_condition(patch, boundary.value, boundary.boundary_condition);
        }
    }
}

// ============================================== Print ============================================
// *************************************************************************************************
/// \brief  Prints info of boundary data
// *************************************************************************************************
void BoundaryDataController::print() const {
#ifndef BENCHMARKING
    for (size_t i = 0; i < number_of_field_types; i++) {
        if (!m_boundary_data[i].is_empty()) {
            m_logger->info("--- found boundary conditions for field {} ({}): ",
                           Mapping::get_field_type_name(FieldType(i)), i);
            m_boundary_data[i].print();
        }
    }
#endif
}

//======================================== Apply boundary condition ================================
// *************************************************************************************************
/// \brief  Applies boundary condition for domain boundary if the field is needed
/// \param  field         Field
/// \param  index_fields  List of indices for each patch
/// \param  patch_starts  List of start indices
/// \param  patch_ends    List of end indices
/// \param  sync          synchronous kernel launching (true, default: false)
// *************************************************************************************************
void BoundaryDataController::apply_boundary_condition(
        Field &field,
        SingleJoinedList **index_fields,
        bool sync) {
    FieldType field_type = field.get_type();
    if (!m_boundary_data[field_type].is_empty()) {
        DomainBoundary::apply_boundary_condition(field,
                                                 index_fields,
                                                 m_boundary_data[field_type],
                                                 sync);
    }
}

//================================ Apply surface boundary condition ================================
// *************************************************************************************************
/// \brief  Applies boundary condition for surfaces in domain boundary if the field is used
/// \param  field         Field
/// \param  index_fields  List of indices for each patch
/// \param  patch_starts  List of start indices
/// \param  patch_ends    List of end indices
/// \param  sync          synchronous kernel launching (true, default: false)
// *************************************************************************************************
void BoundaryDataController::apply_boundary_condition_surface(
        Field &field,
        MultipleJoinedList **index_fields,
        size_t id,
        bool sync) {
    FieldType field_type = field.get_type();
    if (!m_boundary_data[field_type].is_empty()) {
        SurfaceBoundary::apply_boundary_condition(
                field,
                index_fields,
                m_boundary_data[field_type],
                id,
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
        Field &field,
        MultipleJoinedList **index_fields,
        size_t id,
        bool sync) {
    FieldType field_type = field.get_type();
    if (!m_boundary_data[field_type].is_empty()) {
#ifndef BENCHMARKING
        m_logger->debug("apply obstacle boundary conditions of id={} {}", id,
                        Mapping::get_field_type_name(field_type));
#endif
        ObstacleBoundary::apply_boundary_condition(field,
                                                   index_fields,
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
    v_fields.reserve(number_of_field_types);
    for (size_t field_type = 0; field_type < number_of_field_types; field_type++) {
        if (!m_boundary_data[field_type].is_empty()) {
            v_fields.push_back(FieldType(field_type));
        }
    }
    v_fields.shrink_to_fit();
    return v_fields;
}
