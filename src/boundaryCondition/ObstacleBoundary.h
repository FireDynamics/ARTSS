/// \file       ObstacleBoundary.h
/// \brief      Applies boundary condition for obstacle boundary
/// \date       Feb 03, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARYCONDITION_OBSTACLEBOUNDARY_H_
#define ARTSS_BOUNDARYCONDITION_OBSTACLEBOUNDARY_H_

#include <cstdlib>
#include "../utility/GlobalMacrosTypes.h"
#include "../boundary/BoundaryData.h"

namespace ObstacleBoundary {
    void apply_boundary_condition(Settings const &settings,
                                  real* data, size_t** index_fields, const size_t* patch_starts,
                                  const size_t* patch_ends, size_t level,
                                  BoundaryData* boundary_data, size_t id, bool sync = true);
}  // namespace ObstacleBoundary

#endif /* ARTSS_BOUNDARYCONDITION_OBSTACLEBOUNDARY_H_ */

