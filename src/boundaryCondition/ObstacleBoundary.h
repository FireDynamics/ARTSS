/// \file       ObstacleBoundary.h
/// \brief      Applies boundary condition for obstacle boundary
/// \date       Feb 03, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARYCONDITION_OBSTACLEBOUNDARY_H
#define ARTSS_BOUNDARYCONDITION_OBSTACLEBOUNDARY_H

#include <cstdlib>
#include "../utility/GlobalMacrosTypes.h"
#include "../boundary/BoundaryData.h"

namespace ObstacleBoundary {
    void apply_boundary_condition(real* data, size_t** index_fields, const size_t* patch_starts, const size_t* patch_ends, size_t level, BoundaryData* boundary_data, size_t id, bool sync = true);
    void apply_boundary_condition(real* data_field, const size_t* d_patch, size_t patch_start, size_t patch_end, size_t level, int8_t sign_reference_index, size_t reference_index, real value, int8_t sign);
    void apply_neumann(real* data_field, size_t* d_patch, Patch p, size_t patch_start, size_t patch_end, size_t level, real value);
    void apply_dirichlet(real* data_field, size_t* d_patch, Patch p, size_t patch_start, size_t patch_end, size_t level, real value);
    void apply_periodic(real* data_field, size_t* d_patch, Patch p, size_t patch_start, size_t patch_end, size_t level, size_t id);
}  // namespace ObstacleBoundary

#endif /* ARTSS_BOUNDARYCONDITION_OBSTACLEBOUNDARY_H */
