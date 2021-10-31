/// \file       ObstacleBoundary.cpp
/// \brief      Applies boundary condition for obstacle boundary
/// \date       Feb 03, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "ObstacleBoundary.h"
#include "../DomainData.h"
#include "../boundary/BoundaryController.h"


namespace ObstacleBoundary {
namespace {
    //======================================== Apply boundary condition ============================
    // *********************************************************************************************
    /// \brief  Set boundary condition for obstacle boundary
    /// \param  data_field   Field
    /// \param  d_patch List of indices for given patch
    /// \param  patch_start Start Index of Patch
    /// \param  patch_end End index of patch
    /// \param  sign_reference_index Sign of reference index (POSITIVE_SIGN or NEGATIVE_SIGN)
    /// \param  reference_index Index of reference
    /// \param  value Value of boundary condition
    /// \param  sign Sign of boundary condition (POSITIVE_SIGN or NEGATIVE_SIGN)
    // *********************************************************************************************
    void apply_boundary_condition(Field &field, MultipleJoinedList *mjl, size_t id,
                                  int8_t sign_reference_index, size_t reference_index,
                                  real value, int8_t sign) {
        size_t *d_patch = mjl->get_data();
        size_t patch_start = mjl->get_first_index(field.get_level(), id);
        size_t patch_end = mjl->get_last_index(field.get_level(), id);
#ifdef GPU_DEBUG
        auto gpu_logger = Utility::create_gpu_logger("ObstacleBoundary_GPU");
        gpu_logger->debug("applying for [{};{}) with length {} at level {}, pointer {}, field pointer {}",
                      patch_start, patch_end,patch_end-patch_start, field.get_level(),
                      static_cast<void *>(field.data), static_cast<void *>(&field));
        gpu_logger->debug("patch pointer: {}", static_cast<void *>(d_patch));
#endif
#pragma acc data present(field)
        {
#pragma acc parallel loop independent present(d_patch[patch_start:(patch_end-patch_start)]) async
            for (size_t j = patch_start; j < patch_end; ++j) {
                const size_t index = d_patch[j];
                field[index] = sign
                        * field[index + sign_reference_index * static_cast<int>(reference_index)]
                        + value;
            }
#pragma acc wait
        }
    }

    //======================================== Apply dirichlet =====================================
    // *********************************************************************************************
    /// \brief  Apply dirichlet boundary condition
    /// \param  data_field   Field
    /// \param  d_patch List of indices for given patch
    /// \param  p Patch
    /// \param  patch_start Start Index of Patch
    /// \param  patch_end End index of patch
    /// \param  level Multigrid level
    /// \param  value Value of boundary condition
    // *********************************************************************************************
    void apply_dirichlet(Field &field, MultipleJoinedList *mjl, size_t id, Patch p, real value) {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger("ObstacleBoundary");
        logger->debug("applying dirichlet to id={} patch {}", id,
                      PatchObject::get_patch_name(static_cast<Patch>(p)));
#endif
        size_t level = field.get_level();
        if (level > 0) {
            value = 0;
        }
        DomainData *domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx(level);
        size_t Ny = domain->get_Ny(level);
        size_t reference_index = 0;
        int8_t sign_reference_index = POSITIVE_SIGN;
        switch (p) {
            case FRONT:
            case BACK:
                reference_index = Nx * Ny;
                break;

            case BOTTOM:
            case TOP:
                reference_index = Nx;
                break;

            case LEFT:
            case RIGHT:
                reference_index = 1;
                break;

            default:
#ifndef BENCHMARKING
                logger->error("Unknown Patch for dirichlet boundary condition: {}", p);
#endif
                break;
        }

        if (p == FRONT || p == BOTTOM || p == LEFT) {
            sign_reference_index = NEGATIVE_SIGN;
        }
        apply_boundary_condition(field, mjl, id,
                                 sign_reference_index, reference_index, value * 2,
                                 NEGATIVE_SIGN);
    }

    //======================================== apply neumann =======================================
    // *********************************************************************************************
    /// \brief  Apply neumann boundary condition
    /// \param  data_field   Field
    /// \param  d_patch List of indices for given patch
    /// \param  p Patch
    /// \param  patch_start Start Index of Patch
    /// \param  patch_end End index of patch
    /// \param  level Multigrid level
    /// \param  value Value of boundary condition
    // *********************************************************************************************
    void apply_neumann(Field &field, MultipleJoinedList *mjl, size_t id, Patch p, real value) {
        size_t level = field.get_level();
        if (level > 0) {
            value = 0;
        }
        DomainData *domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx(level);
        size_t Ny = domain->get_Ny(level);
        size_t reference_index = 0;
        int8_t sign_reference_index = POSITIVE_SIGN;
        switch (p) {
            case FRONT:
            case BACK:
                value *= domain->get_dz(level);
                reference_index = Nx * Ny;
                break;
            case TOP:
            case BOTTOM:
                value *= domain->get_dy(level);
                reference_index = Nx;
                break;
            case LEFT:
            case RIGHT:
                value *= domain->get_dx(level);
                reference_index = 1;
                break;
            default:
#ifndef BENCHMARKING
                auto logger = Utility::create_logger("ObstacleBoundary");
                logger->error("Unknown Patch for neumann boundary condition: {}", p);
#endif
                break;
        }

        if (p == FRONT || p == BOTTOM || p == LEFT) {
            sign_reference_index = NEGATIVE_SIGN;
        }

        apply_boundary_condition(field, mjl, id,
                                 sign_reference_index, reference_index, -value, POSITIVE_SIGN);
    }
}  // namespace

//======================================== Apply boundary condition ================================
// *************************************************************************************************
/// \brief  Applies boundary condition for obstacle boundary
/// \param  field   Field
/// \param  index_fields List of indices for each patch
/// \param  patch_starts List of start indices
/// \param  patch_ends List of end indices
/// \param  boundary_data Boundary data object of Domain
/// \param  id ID of obstacle
/// \param  sync synchronous kernel launching (true, default: false)
// *************************************************************************************************
void apply_boundary_condition(Field &field, MultipleJoinedList **index_fields,
                              BoundaryData *boundary_data, size_t id, bool sync) {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger("ObstacleBoundary");
#endif
        size_t level = field.get_level();
        for (size_t p = 0; p < number_of_patches; p++) {
            auto patch = static_cast<Patch>(p);
            size_t size_patch = index_fields[patch]->get_slice_size(level, id);
            if (size_patch == 0) {
#ifndef BENCHMARKING
                logger->debug("skipping apply boundary condition id={} for: {}",
                              id,PatchObject::get_patch_name(patch));
#endif
                continue;
            }
            BoundaryCondition bc = boundary_data->get_boundary_condition(patch);
            switch (bc) {
                case BoundaryCondition::DIRICHLET:
                    apply_dirichlet(field, index_fields[patch], id, patch,
                                    boundary_data->get_value(patch));
                    break;
                case BoundaryCondition::NEUMANN:
                    apply_neumann(field, index_fields[patch], id, patch,
                                  boundary_data->get_value(patch));
                    break;
                case BoundaryCondition::PERIODIC:
#ifndef BENCHMARKING
                    logger->error("periodic boundary conditions are not implemented for obstacles");
#endif
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

}  // namespace ObstacleBoundary
