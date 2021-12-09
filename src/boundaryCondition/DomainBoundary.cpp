/// \file       DomainBoundary.cpp
/// \brief      Applies boundary condition for domain boundary
/// \date       Feb 03, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "DomainBoundary.h"
#include "../DomainData.h"

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
        void apply_boundary_condition(Settings::Settings const &settings,
                Field &field, SingleJoinedList *jl,
                int8_t sign_reference_index, size_t reference_index, real value, int8_t sign) {
            size_t level = field.get_level();
            real *data = field.data;
            size_t patch_start = jl->get_first_index(level);
            size_t patch_end = jl->get_last_index(level);
            size_t patch_size = jl->get_slice_size(level);
            size_t *d_patch = jl->get_data();
#ifndef BENCHMARKING
            auto logger = Utility::create_logger(settings, "DomainBoundary");
            logger->debug("apply_boundary_condition ! apply domain boundary pointer patch: {} size: {}", static_cast<void *>(d_patch), jl->get_size());
            logger->debug("apply_boundary_condition ! apply domain boundary pointer field: {} size: {}", static_cast<void *>(field.data), field.get_size());
            logger->debug("apply_boundary_condition ! start {} end {} size {} level {}", patch_start, patch_end, jl->get_slice_size(level), level);
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
        /// \param  jl single joined list which contains the patch data
        /// \param  patch Patch
        /// \param  value Value of boundary condition
        // *********************************************************************************************
        void apply_dirichlet(Settings::Settings const &settings, Field &field, SingleJoinedList *jl, Patch patch, real value) {
            DomainData *domain = DomainData::getInstance();
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

            if (patch == BACK || patch == TOP || patch == RIGHT) {
                sign_reference_index = NEGATIVE_SIGN;
            }

            apply_boundary_condition(settings, field, jl, sign_reference_index, reference_index, value * 2, NEGATIVE_SIGN);
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
        void apply_neumann(Settings::Settings const &settings, Field &field, SingleJoinedList *jl, Patch patch, real value) {
            size_t level = field.get_level();
            DomainData *domain = DomainData::getInstance();
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

            if (patch == BACK || patch == TOP || patch == RIGHT) {
                sign_reference_index = NEGATIVE_SIGN;
            }
            apply_boundary_condition(settings, field, jl, sign_reference_index, reference_index, value, POSITIVE_SIGN);
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
        void apply_periodic(Settings::Settings const &settings, Field &field, SingleJoinedList *jl, Patch patch) {
            size_t level = field.get_level();
            DomainData *domain = DomainData::getInstance();
            size_t Nx = domain->get_Nx(level);
            size_t Ny = domain->get_Ny(level);

            size_t reference_index = 0;
            int8_t sign_reference_index = POSITIVE_SIGN;

            switch (patch) {
                case FRONT:
                case BACK:
                    reference_index = Nx * Ny * DomainData::getInstance()->get_nz(level);
                    break;
                case BOTTOM:
                case TOP:
                    reference_index = Nx * DomainData::getInstance()->get_ny(level);
                    break;
                case LEFT:
                case RIGHT:
                    reference_index = DomainData::getInstance()->get_nx(level);
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

            apply_boundary_condition(settings, field, jl, sign_reference_index, reference_index, 0, POSITIVE_SIGN);
        }
    }

//======================================== Apply boundary condition ================================
// *************************************************************************************************
/// \brief  Applies boundary condition for domain boundary
/// \param  data_field   Field
/// \param  index_fields List of indices for each patch
/// \param  boundary_data Boundary data_field object of Domain
/// \param  sync synchronous kernel launching (true, default: false)
// *************************************************************************************************
    void apply_boundary_condition(Settings::Settings const &settings, Field &field, SingleJoinedList **index_fields, BoundaryData *boundary_data, bool sync) {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(settings, "DomainBoundary");
#endif
        for (size_t i = 0; i < number_of_patches; i++) {
            size_t level = field.get_level();
            SingleJoinedList *jl = index_fields[i];
            size_t patch_size = jl->get_slice_size(level);
            if (patch_size == 0) {
                continue;
            }
            auto p = static_cast<Patch>(i);
#ifndef BENCHMARKING
            logger->debug("apply_boundary_condition ! level {} for {}", jl->get_slice_size(level), level, PatchObject::get_patch_name(p));
#endif
            BoundaryCondition bc = boundary_data->get_boundary_condition(p);
            real value = 0;
            switch (bc) {
                case BoundaryCondition::DIRICHLET:
                    if (field.get_level() == 0) {
                        value = boundary_data->get_value(p);
                    }
                    apply_dirichlet(settings, field, jl, p, value);
                    break;
                case BoundaryCondition::NEUMANN:
                    if (field.get_level() == 0) {
                        value = boundary_data->get_value(p);
                    }
                    apply_neumann(settings, field, jl, p, value);
                    break;
                case BoundaryCondition::PERIODIC:
                    apply_periodic(settings, field, jl, p);
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
}  // namespace DomainBoundary
