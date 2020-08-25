/// \file       ObstacleBoundary.cpp
/// \brief      Applies boundary condition for obstacle boundary
/// \date       Feb 03, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "ObstacleBoundary.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"

//======================================== Apply boundary condition ====================================
// ***************************************************************************************
/// \brief  Applies boundary condition for obstacle boundary
/// \param  dataField   Field
/// \param  index_fields List of indices for each patch
/// \param  patch_starts List of start indices
/// \param  patch_ends List of end indices
/// \param  level Multigrid level
/// \param  boundary_data Boundary data object of Domain
/// \param  id ID of obstacle
/// \param  sync synchronous kernel launching (true, default: false)
// ***************************************************************************************
void ObstacleBoundary::applyBoundaryCondition(real *dataField, size_t **index_fields, const size_t *patch_starts, const size_t *patch_ends, size_t level, BoundaryData *boundary_data, size_t id, bool sync) {
    for (size_t i = 0; i < numberOfPatches; i++) {
        size_t *d_patch = *(index_fields + i);
        size_t patch_start = *(patch_starts + i);
        size_t patch_end = *(patch_ends + i);
        Patch p = static_cast<Patch>(i);
        BoundaryCondition bc = boundary_data->getBoundaryCondition(p);
        switch (bc) {
            case BoundaryCondition::DIRICHLET:
                applyDirichlet(dataField, d_patch, p, patch_start, patch_end, level, boundary_data->getValue(p));
                break;
            case BoundaryCondition::NEUMANN:
                applyNeumann(dataField, d_patch, p, patch_start, patch_end, level, boundary_data->getValue(p));
                break;
            case BoundaryCondition::PERIODIC:
                applyPeriodic(dataField, d_patch, p, patch_start, patch_end, level, id);
                break;
            default:
#ifndef BENCHMARKING
                auto m_logger = Utility::create_logger("ObstacleBoundary");
                m_logger->error("Unknown boundary condition: {}", bc);
#endif
                break;
        }
    }
    if (sync) {
#pragma acc wait
    }
}

//======================================== Apply boundary condition ====================================
// ***************************************************************************************
/// \brief  Set boundary condition for obstacle boundary
/// \param  data_field   Field
/// \param  d_patch List of indices for given patch
/// \param  patch_start Start Index of Patch
/// \param  patch_end End index of patch
/// \param  level Multigrid level
/// \param  sign_reference_index Sign of reference index (POSITIVE_SIGN or NEGATIVE_SIGN)
/// \param  reference_index Index of reference
/// \param  value Value of boundary condition
/// \param  sign Sign of boundary condition (POSITIVE_SIGN or NEGATIVE_SIGN)
// ***************************************************************************************
void ObstacleBoundary::applyBoundaryCondition(real *data_field, const size_t *d_patch, size_t patch_start, size_t patch_end, size_t level, int8_t sign_reference_index, size_t reference_index, real value, int8_t sign) {
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
void ObstacleBoundary::applyDirichlet(real *data_field, size_t *d_patch, Patch patch, size_t patch_start, size_t patch_end, size_t level, real value) {
    if (level > 0) {
        value = 0;
    }
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);
    size_t reference_index = 0;
    int8_t sign_reference_index = POSITIVE_SIGN;
    switch (patch) {
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
            auto m_logger = Utility::create_logger("ObstacleBoundary");
            m_logger->error("Unknown Patch for dirichlet boundary condition: {}", patch);
#endif
            break;
    }

    if (patch == FRONT || patch == BOTTOM || patch == LEFT){
        sign_reference_index = NEGATIVE_SIGN;
    }
    applyBoundaryCondition(
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
void ObstacleBoundary::applyNeumann(real *data_field, size_t *d_patch, Patch patch, size_t patch_start, size_t patch_end, size_t level, real value) {
    if (level > 0) {
        value = 0;
    }
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);
    size_t reference_index = 0;
    int8_t sign_reference_index = POSITIVE_SIGN;
    switch (patch) {
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
            auto m_logger = Utility::create_logger("ObstacleBoundary");
            m_logger->error("Unknown Patch for neumann boundary condition: {}", patch);
#endif
            break;
    }

    if (patch == FRONT || patch == BOTTOM || patch == LEFT){
        sign_reference_index = NEGATIVE_SIGN;
    }

    applyBoundaryCondition(
            data_field, d_patch, patch_start, patch_end,
            level, sign_reference_index, reference_index, -value, POSITIVE_SIGN);
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
void ObstacleBoundary::applyPeriodic(real *data_field, size_t *d_patch, Patch patch, size_t patch_start, size_t patch_end, size_t level, size_t id) {
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);
    BoundaryController *bdc = BoundaryController::getInstance();

    size_t reference_index = 0;
    int8_t sign_reference_index = POSITIVE_SIGN;
    switch (patch) {
        case FRONT:
        case BACK:
            reference_index = Nx * Ny * (bdc->getObstacleStrideZ(id, level) - 2);
            break;
        case BOTTOM:
        case TOP:
            reference_index = Nx * (bdc->getObstacleStrideY(id, level) - 2);
            break;
        case LEFT:
        case RIGHT:
            reference_index = (bdc->getObstacleStrideX(id, level) - 2);
            break;
        default:
#ifndef BENCHMARKING
            auto m_logger = Utility::create_logger("ObstacleBoundary");
            m_logger->error("Unknown Patch for periodic boundary condition: {}", patch);
#endif
            break;
    }

    if (patch == BACK || patch == TOP || patch == RIGHT){
        sign_reference_index = NEGATIVE_SIGN;
    }
    applyBoundaryCondition(
            data_field, d_patch, patch_start, patch_end,
            level, sign_reference_index, reference_index, 0, POSITIVE_SIGN);
}
