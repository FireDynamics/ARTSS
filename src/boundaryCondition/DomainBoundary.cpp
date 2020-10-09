/// \file       DomainBoundary.cpp
/// \brief      Applies boundary condition for domain boundary
/// \date       Feb 03, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "DomainBoundary.h"
#include "../Domain.h"

//======================================== Apply boundary condition ====================================
// ***************************************************************************************
/// \brief  Applies boundary condition for domain boundary
/// \param  data_field   Field
/// \param  index_fields List of indices for each patch
/// \param  patch_start List of start indices
/// \param  patch_end List of end indices
/// \param  level Multigrid level
/// \param  boundary_data Boundary data_field object of Domain
/// \param  sync synchronous kernel launching (true, default: false)
// ***************************************************************************************
void DomainBoundary::apply_boundary_condition(real *data_field, size_t **index_fields, const size_t *patch_starts, const size_t *patch_ends, size_t level, BoundaryData *boundary_data, bool sync) {

#ifdef USEMPI
    auto mpi_handler = MPIHandler::getInstance();
    std::vector<int> rank_has_neighbour{mpi_handler->get_mpi_neighbour()};
#endif

    for (size_t i = 0; i < numberOfPatches; i++) {
        size_t *d_patch = *(index_fields + i);
        size_t patch_start = *(patch_starts + i);
        size_t patch_end = *(patch_ends + i);
        Patch p = static_cast<Patch>(i);
        BoundaryCondition bc = boundary_data->getBoundaryCondition(p);
        switch (bc) {
            case BoundaryCondition::DIRICHLET:
                apply_dirichlet(data_field, d_patch, p, patch_start, patch_end, level, boundary_data->getValue(p));
                break;
            case BoundaryCondition::NEUMANN:
                apply_neumann(data_field, d_patch, p, patch_start, patch_end, level, boundary_data->getValue(p));
                break;
            case BoundaryCondition::PERIODIC:
                apply_periodic(data_field, d_patch, p, patch_start, patch_end, level);
                break;
            default:
#ifndef BENCHMARKING
                auto logger = Utility::create_logger("DomainBoundary");
                logger->error("Unknown boundary condition: {}", bc);
#endif
                break;
        }
#ifdef USEMPI
        if(rank_has_neighbour.at(i) == 1) {
            mpi_handler->exchange_data(data_field, p, d_patch, patch_start, level);
            mpi_handler->set_barrier();
        }
#endif
    }

    if (sync) {
#pragma acc wait
    }
}

//======================================== Apply boundary condition ====================================
// ***************************************************************************************
/// \brief  Set boundary condition for domain boundary
/// \param  data_field   Field
/// \param  d_patch List of indices for given patch
/// \param  patch_start Start Index of Patch
/// \param  patch_end End index of patch
/// \param  level Multigrid level
/// \param  sign_reference_index Sign of reference index ( POSITIVE_SIGN or NEGATIVE_SIGN )
/// \param  reference_index Index of reference
/// \param  value Value of boundary condition
/// \param  sign Sign of boundary condition ( POSITIVE_SIGN or NEGATIVE_SIGN )
// ***************************************************************************************
void DomainBoundary::apply_boundary_condition(real *data_field, const size_t *d_patch, size_t patch_start, size_t patch_end, size_t level, int8_t sign_reference_index, size_t reference_index, real value, int8_t sign) {
    Domain *domain = Domain::getInstance();
    size_t b_size = domain->get_size(level);
#pragma acc data present(data_field[:b_size])
    {
#pragma acc parallel loop independent present(d_patch[patch_start:(patch_end-patch_start)]) async
        for (size_t j = patch_start; j < patch_end; ++j) {
            const size_t index = d_patch[j];
            data_field[index] = sign * data_field[index + sign_reference_index * static_cast<int>(reference_index)] + value;
        }
#pragma acc wait
    }
}

//======================================== Apply dirichlet ====================================
// ***************************************************************************************
/// \brief  Apply dirichlet boundary condition
/// \param  data_field   Field
/// \param  d_patch List of indices for given patch
/// \param  patch Patch
/// \param  patch_start Start Index of Patch
/// \param  patch_end End index of patch
/// \param  level Multigrid level
/// \param  value Value of boundary condition
// ***************************************************************************************
void DomainBoundary::apply_dirichlet(real *data_field, size_t *d_patch, Patch patch, size_t patch_start, size_t patch_end, size_t level, real value) {
    if (level > 0) {
        value = 0;
    }

    Domain *domain = Domain::getInstance();
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
            auto logger = Utility::create_logger("DomainBoundary");
            logger->error("Unknown Patch for dirichlet boundary condition: {}", patch);
#endif
            break;
    }

    if (patch == BACK || patch == TOP || patch == RIGHT){
        sign_reference_index = NEGATIVE_SIGN;
    }

    apply_boundary_condition(
            data_field, d_patch, patch_start, patch_end,
            level, sign_reference_index, reference_index, value * 2, NEGATIVE_SIGN);
}

//======================================== Apply neumann ====================================
// ***************************************************************************************
/// \brief  Apply neumann boundary condition
/// \param  data_field   Field
/// \param  d_patch List of indices for given patch
/// \param  patch Patch
/// \param  patch_start Start Index of Patch
/// \param  patch_end End index of patch
/// \param  level Multigrid level
/// \param  value Value of boundary condition
// ***************************************************************************************
void DomainBoundary::apply_neumann(real *data_field, size_t *d_patch, Patch patch, size_t patch_start, size_t patch_end, size_t level, real value) {
    if (level > 0) {
        value = 0;
    }
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
            auto logger = Utility::create_logger("DomainBoundary");
            logger->error("Unknown Patch for neumann boundary condition: {}", patch);
#endif
            break;
    }

    if (patch == BACK || patch == TOP || patch == RIGHT){
        sign_reference_index = NEGATIVE_SIGN;
    }
    apply_boundary_condition(
            data_field, d_patch, patch_start, patch_end,
            level, sign_reference_index, reference_index, value, POSITIVE_SIGN);
}

//======================================== Apply periodic ====================================
// ***************************************************************************************
/// \brief  Apply periodic boundary condition
/// \param  data_field   Field
/// \param  d_patch List of indices for given patch
/// \param  patch Patch
/// \param  patch_start Start Index of Patch
/// \param  patch_end End index of patch
/// \param  level Multigrid level
// ***************************************************************************************
void DomainBoundary::apply_periodic(real *data_field, size_t *d_patch, Patch patch, size_t patch_start, size_t patch_end, size_t level) {
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);

    size_t reference_index = 0;
    int sign_reference_index = POSITIVE_SIGN;

    switch (patch) {
        case FRONT:
        case BACK:
            reference_index = Nx * Ny * (Domain::getInstance()->get_nz(level) - 2);
            break;
        case BOTTOM:
        case TOP:
            reference_index = Nx * (Domain::getInstance()->get_ny(level) - 2);
            break;
        case LEFT:
        case RIGHT:
            reference_index = (Domain::getInstance()->get_nx(level) - 2);
            break;
        default:
#ifndef BENCHMARKING
            auto logger = Utility::create_logger("DomainBoundary");
            logger->error("Unknown Patch for periodic boundary condition: {}", patch);
#endif
            break;
    }

    if (patch == BACK || patch == TOP || patch == RIGHT) {
        sign_reference_index = NEGATIVE_SIGN;
    }

    apply_boundary_condition(
            data_field, d_patch, patch_start, patch_end,
            level, sign_reference_index, reference_index, 0, POSITIVE_SIGN);
}
