/// \file       SurfaceBoundary.h
/// \brief      Applies boundary condition for surface boundary
/// \date       Dec 01, 2021
/// \author     My Linh Würzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARYCONDITION_SURFACEBOUNDARY_H
#define ARTSS_BOUNDARYCONDITION_SURFACEBOUNDARY_H

#include "../domain/BoundaryData.h"
#include "../field/Field.h"
#include "../GPULists/MultipleJoinedList.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"

namespace SurfaceBoundary {
    void apply_boundary_condition(Field &field, MultipleJoinedList** index_fields,
                                  const BoundaryData &boundary_data, size_t id, bool sync = true);
}  // namespace SurfaceBoundary
#endif /* ARTSS_BOUNDARYCONDITION_SURFACEBOUNDARY_H */
