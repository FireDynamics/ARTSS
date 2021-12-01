/// \file       DomainBoundary.cpp
/// \brief      Applies boundary condition for domain boundary
/// \date       Feb 03, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "DomainBoundary.h"
#include "../Domain.h"

namespace DomainBoundary {
namespace {
    //======================================== Apply boundary condition ============================
    // *********************************************************************************************
    /// \brief  Set boundary condition for domain boundary
    /// \param  field   Field
    /// \param  d_patch List of indices for given patch
    /// \param  patch_start Start Index of Patch
    /// \param  patch_end End index of patch
    /// \param  sign_reference_index Sign of reference index ( POSITIVE_SIGN or NEGATIVE_SIGN )
    /// \param  reference_index Index of reference
    /// \param  value Value of boundary condition
    /// \param  sign Sign of boundary condition ( POSITIVE_SIGN or NEGATIVE_SIGN )
    // *********************************************************************************************
    void apply_boundary_condition(
            Field &field, const size_t *d_patch, const size_t patch_start, const size_t patch_end,
            int8_t sign_reference_index, size_t reference_index, real value, int8_t sign) {
        real *data = field.data;
#pragma acc data present(field)
        {
#pragma acc parallel loop independent present(d_patch[patch_start:(patch_end-patch_start)]) async
            for (size_t j = patch_start; j < patch_end; ++j) {
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
    /// \param  d_patch List of indices for given patch
    /// \param  patch Patch
    /// \param  patch_start Start Index of Patch
    /// \param  patch_end End index of patch
    /// \param  value Value of boundary condition
    // *********************************************************************************************
    void apply_dirichlet(Settings::Settings const &settings, Field &field, size_t *d_patch, Patch patch,
                         const size_t patch_start, const size_t patch_end, real value) {
        Domain *domain = Domain::getInstance();
        size_t level = field.get_level();
        size_t reference_index = 0;
        int8_t sign_reference_index = POSITIVE_SIGN;
        switch (patch) {
            case FRONT:
            case BACK:
                reference_index = domain->get_Nx(level) * domain->get_Ny(level);
                break;
            case BOTTOM:
            case TOP:
                reference_index = domain->get_Nx(level);
                break;
            case LEFT:
            case RIGHT:
                reference_index = 1;
                break;
            default:
#ifndef BENCHMARKING
                auto logger = Utility::create_logger(settings, "DomainBoundary");
                logger->error("Unknown Patch for dirichlet boundary condition: {}", patch);
#endif
                break;
        }

        if (patch == BACK || patch == TOP || patch == RIGHT){
            sign_reference_index = NEGATIVE_SIGN;
        }

        apply_boundary_condition(field, d_patch, patch_start, patch_end,
                                 sign_reference_index, reference_index, value * 2, NEGATIVE_SIGN);
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
    void apply_neumann(Settings::Settings const &settings, Field &field, size_t *d_patch, Patch patch,
                       size_t patch_start, size_t patch_end, real value) {
        size_t level = field.get_level();
        Domain *domain = Domain::getInstance();
        size_t reference_index = 0;
        int8_t sign_reference_index = POSITIVE_SIGN;
        switch (patch) {
            case FRONT:
            case BACK:
                value *= domain->get_dz(level);
                reference_index = domain->get_Nx(level) * domain->get_Ny(level);
                break;
            case BOTTOM:
            case TOP:
                value *= domain->get_dy(level);
                reference_index = domain->get_Nx(level);
                break;
            case LEFT:
            case RIGHT:
                value *= domain->get_dz(level);
                reference_index = 1;
                break;
            default:
#ifndef BENCHMARKING
                auto logger = Utility::create_logger(settings, "DomainBoundary");
                logger->error("Unknown Patch for neumann boundary condition: {}", patch);
#endif
                break;
        }

        if (patch == BACK || patch == TOP || patch == RIGHT){
            sign_reference_index = NEGATIVE_SIGN;
        }
        apply_boundary_condition(field, d_patch, patch_start, patch_end,
                                 sign_reference_index, reference_index, value, POSITIVE_SIGN);
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
    void apply_periodic(Settings::Settings const &settings, Field &field, size_t *d_patch, Patch patch,
                        const size_t patch_start, const size_t patch_end) {
        size_t level = field.get_level();
        Domain *domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(level);
        size_t Ny = domain->get_Ny(level);

        size_t reference_index = 0;
        int8_t sign_reference_index = POSITIVE_SIGN;

        switch (patch) {
            case FRONT:
            case BACK:
                reference_index = Nx * Ny * Domain::getInstance()->get_nz(level);
                break;
            case BOTTOM:
            case TOP:
                reference_index = Nx * Domain::getInstance()->get_ny(level);
                break;
            case LEFT:
            case RIGHT:
                reference_index = Domain::getInstance()->get_nx(level);
                break;
            default:
#ifndef BENCHMARKING
                auto logger = Utility::create_logger(settings, "DomainBoundary");
                logger->error("Unknown Patch for periodic boundary condition: {}", patch);
#endif
                break;
        }

        if (patch == BACK || patch == TOP || patch == RIGHT) {
            sign_reference_index = NEGATIVE_SIGN;
        }

        apply_boundary_condition(field, d_patch, patch_start, patch_end,
                                 sign_reference_index, reference_index, 0, POSITIVE_SIGN);
    }
}

//======================================== Apply boundary condition ================================
// *************************************************************************************************
/// \brief  Applies boundary condition for domain boundary
/// \param  data_field   Field
/// \param  index_fields List of indices for each patch
/// \param  patch_start List of start indices
/// \param  patch_end List of end indices
/// \param  level Multigrid level
/// \param  boundary_data Boundary data_field object of Domain
/// \param  sync synchronous kernel launching (true, default: false)
// *************************************************************************************************
void apply_boundary_condition(Settings::Settings const &settings, Field &field, size_t **index_fields, const size_t *patch_starts,
                              const size_t *patch_ends, BoundaryData *boundary_data, bool sync) {
    for (size_t i = 0; i < number_of_patches; i++) {
        size_t *d_patch = *(index_fields + i);
        size_t patch_start = *(patch_starts + i);
        size_t patch_end = *(patch_ends + i);
        auto p = static_cast<Patch>(i);
        BoundaryCondition bc = boundary_data->get_boundary_condition(p);
        real value = 0;
        switch (bc) {
            case BoundaryCondition::DIRICHLET:
                if (field.get_level() == 0) {
                    value = boundary_data->get_value(p);
                }
                apply_dirichlet(settings, field, d_patch, p, patch_start, patch_end, value);
                break;
            case BoundaryCondition::NEUMANN:
                if (field.get_level() == 0) {
                    value = boundary_data->get_value(p);
                }
                apply_neumann(settings, field, d_patch, p, patch_start, patch_end, value);
                break;
            case BoundaryCondition::PERIODIC:
                apply_periodic(settings, field, d_patch, p, patch_start, patch_end);
                break;
            default:
#ifndef BENCHMARKING
                auto logger = Utility::create_logger(settings, "DomainBoundary");
                logger->error("Unknown boundary condition: {}", bc);
#endif
                break;
        }
    }
    if (sync) {
#pragma acc wait
    }
}
}  // namespace DomainBoundary
