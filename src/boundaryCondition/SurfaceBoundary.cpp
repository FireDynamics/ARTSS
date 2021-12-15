/// \file       SurfaceBoundary.cpp
/// \brief      Applies boundary condition for surface boundary
/// \date       Dec 01, 2021
/// \author     My Linh Würzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.

#include "SurfaceBoundary.h"
#include "../boundary/DomainData.h"

namespace SurfaceBoundary {
    std::string class_name = "SurfaceBoundary";
namespace {
//======================================== Apply boundary condition ============================
// *********************************************************************************************
/// \brief  Set boundary condition for surface boundary
/// \param  field   Field
/// \param  d_patch List of indices for given patch
/// \param  patch_start Start Index of Patch
/// \param  patch_end End index of patch
/// \param  sign_reference_index Sign of reference index ( POSITIVE_SIGN or NEGATIVE_SIGN )
/// \param  reference_index Index of reference
/// \param  value Value of boundary condition
/// \param  sign Sign of boundary condition ( POSITIVE_SIGN or NEGATIVE_SIGN )
// *********************************************************************************************
void apply_boundary_condition(Field &field, MultipleJoinedList *mjl, size_t id,
                              int8_t sign_reference_index, size_t reference_index,
                              real value, int8_t sign) {
    size_t level = field.get_level();
    real *data = field.data;
    size_t patch_start = mjl->get_first_index(level, id);
    size_t patch_end = mjl->get_last_index(level, id);
    size_t patch_size __attribute__((unused)) = mjl->get_slice_size(level, id);
    size_t *d_patch = mjl->get_data();
#ifndef BENCHMARKING
    auto logger = Utility::create_logger(class_name);
    logger->debug("apply_boundary_condition ! apply surface boundary pointer patch: {} size: {}", static_cast<void *>(d_patch), mjl->get_size());
    logger->debug("apply_boundary_condition ! apply surface boundary pointer field: {} size: {}", static_cast<void *>(field.data), field.get_size());
    logger->debug("apply_boundary_condition ! start {} end {} size {} level {}", patch_start, patch_end, mjl->get_slice_size(level), level);
#endif
#pragma acc data present(field)
{
#pragma acc parallel loop independent present(d_patch[patch_start:patch_size]) async
    for (size_t j = patch_start; j <= patch_end; ++j) {
        const size_t index = d_patch[j];
        *(data + index) = sign * *(data + index + sign_reference_index * reference_index) + value;
    }
#pragma acc wait
}
}

//======================================== Apply dirichlet =====================================
// *********************************************************************************************
/// \brief  Apply dirichlet boundary condition
/// \param  field   Field
/// \param  mjl multiple joined list which contains the patch data
/// \param  patch Patch
/// \param  value Value of boundary condition
// *********************************************************************************************
void apply_dirichlet(Field &field, MultipleJoinedList *mjl,
                     size_t id, Patch patch, real value) {
    auto domain_data = DomainData::getInstance();
    size_t level = field.get_level();
    size_t reference_index = 0;
    int8_t sign_reference_index = POSITIVE_SIGN;
    switch (patch) {
        case FRONT:
        case BACK:
            reference_index = domain_data->get_Nx(level) * domain_data->get_Ny(level);
            break;
        case BOTTOM:
        case TOP:
            reference_index = domain_data->get_Nx(level);
            break;
        case LEFT:
        case RIGHT:
            reference_index = 1;
            break;
        default:
#ifndef BENCHMARKING
            auto logger = Utility::create_logger(class_name);
            logger->error("Unknown Patch for dirichlet boundary condition: {}", patch);
#endif
            break;
    }

    if (patch == BACK || patch == TOP || patch == RIGHT) {
        sign_reference_index = NEGATIVE_SIGN;
    }

    apply_boundary_condition(field, mjl, id, sign_reference_index, reference_index, value * 2, NEGATIVE_SIGN);
}

//======================================== Apply neumann =======================================
// *********************************************************************************************
/// \brief  Apply neumann boundary condition
/// \param  field   Field
/// \param  d_patch List of indices for given patch
/// \param  patch Patch
/// \param  patch_start Start Index of Patch
/// \param  patch_end End index of patch
/// \param  value Value of boundary condition
// *********************************************************************************************
void apply_neumann(Field &field, MultipleJoinedList *mjl,
                   size_t id, Patch patch, real value) {
    size_t level = field.get_level();
    auto domain_data = DomainData::getInstance();
    size_t reference_index = 0;
    int8_t sign_reference_index = POSITIVE_SIGN;
    switch (patch) {
        case FRONT:
        case BACK:
            value *= domain_data->get_dz(level);
            reference_index = domain_data->get_Nx(level) * domain_data->get_Ny(level);
            break;
        case BOTTOM:
        case TOP:
            value *= domain_data->get_dy(level);
            reference_index = domain_data->get_Nx(level);
            break;
        case LEFT:
        case RIGHT:
            value *= domain_data->get_dz(level);
            reference_index = 1;
            break;
        default:
#ifndef BENCHMARKING
            auto logger = Utility::create_logger(class_name);
            logger->error("Unknown Patch for neumann boundary condition: {}", patch);
#endif
            break;
    }

    if (patch == BACK || patch == TOP || patch == RIGHT) {
        sign_reference_index = NEGATIVE_SIGN;
    }
    apply_boundary_condition(field, mjl, id, sign_reference_index, reference_index, value, POSITIVE_SIGN);
}

//======================================== Apply periodic ======================================
// *********************************************************************************************
/// \brief  Apply periodic boundary condition
/// \param  data_field   Field
/// \param  d_patch List of indices for given patch
/// \param  patch Patch
/// \param  patch_start Start Index of Patch
/// \param  patch_end End index of patch
/// \param  level Multigrid level
// *********************************************************************************************
void apply_periodic(Field &field, MultipleJoinedList *mjl,
                    size_t id, Patch patch) {
    size_t level = field.get_level();
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx(level);
    size_t Ny = domain_data->get_Ny(level);

    size_t reference_index = 0;
    int8_t sign_reference_index = POSITIVE_SIGN;

    switch (patch) {
        case FRONT:
        case BACK:
            reference_index = Nx * Ny * domain_data->get_nz(level);
            break;
        case BOTTOM:
        case TOP:
            reference_index = Nx * domain_data->get_ny(level);
            break;
        case LEFT:
        case RIGHT:
            reference_index = domain_data->get_nx(level);
            break;
        default:
#ifndef BENCHMARKING
            auto logger = Utility::create_logger(class_name);
            logger->error("Unknown Patch for periodic boundary condition: {}", patch);
#endif
                    break;
            }

            if (patch == BACK || patch == TOP || patch == RIGHT) {
                sign_reference_index = NEGATIVE_SIGN;
            }

            apply_boundary_condition(field, mjl, id, sign_reference_index, reference_index, 0, POSITIVE_SIGN);
        }
    }

//======================================== Apply boundary condition ================================
// *************************************************************************************************
/// \brief  Applies boundary condition for surface boundary
/// \param  data_field   Field
/// \param  index_fields List of indices for each patch
/// \param  boundary_data Boundary data_field object of Domain
/// \param  sync synchronous kernel launching (true, default: false)
// *************************************************************************************************
void apply_boundary_condition(Field &field, MultipleJoinedList **index_fields,
                              const BoundaryData &boundary_data, size_t id, bool sync) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger(class_name);
#endif
    for (size_t i = 0; i < number_of_patches; i++) {
        size_t level = field.get_level();
        MultipleJoinedList *mjl = index_fields[i];
        size_t patch_size = mjl->get_slice_size(level);
        if (patch_size == 0) {
            continue;
        }
        auto p = Patch(i);
#ifndef BENCHMARKING
        logger->debug("apply_boundary_condition ! level {} for {}", mjl->get_slice_size(level), level, Mapping::get_patch_name(p));
#endif
        BoundaryCondition bc = boundary_data.get_boundary_condition(p);
        real value = 0;
        switch (bc) {
            case BoundaryCondition::DIRICHLET:
                if (field.get_level() == 0) {
                    value = boundary_data.get_value(p);
                }
                apply_dirichlet(field, mjl, id, p, value);
                break;
            case BoundaryCondition::NEUMANN:
                if (field.get_level() == 0) {
                    value = boundary_data.get_value(p);
                }
                apply_neumann(field, mjl, id, p, value);
                break;
            case BoundaryCondition::PERIODIC:
                apply_periodic(field, mjl, id, p);
                break;
            default:
#ifndef BENCHMARKING
                logger->error("Unknown boundary condition: {}", bc);
#endif
                break;
        }
    }
    if (sync) {
#pragma acc wait
    }
}
}  // namespace SurfaceBoundary
