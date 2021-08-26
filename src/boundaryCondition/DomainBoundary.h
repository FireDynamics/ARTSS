/// \file       DomainBoundary.h
/// \brief      Applies boundary condition for domain boundary
/// \date       Feb 03, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARYCONDITION_DOMAINBOUNDARY_H
#define ARTSS_BOUNDARYCONDITION_DOMAINBOUNDARY_H

#include "../utility/GlobalMacrosTypes.h"
#include "../boundary/BoundaryData.h"
#include "../utility/Utility.h"

namespace DomainBoundary {
    void apply_boundary_condition(real* data_field, size_t** index_fields,
                                  const size_t* patch_starts, const size_t* patch_ends,
                                  size_t level, BoundaryData* boundary_data, bool sync = true);
}  // namespace DomainBoundary
#endif /* ARTSS_BOUNDARYCONDITION_DOMAINBOUNDARY_H */
