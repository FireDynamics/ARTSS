/// \file       DomainBoundary.h
/// \brief      Applies boundary condition for domain boundary
/// \date       Feb 03, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARYCONDITION_DOMAINBOUNDARY_H
#define ARTSS_BOUNDARYCONDITION_DOMAINBOUNDARY_H

#include "../boundary/BoundaryData.h"
#include "../field/Field.h"
#include "../utility/Utility.h"
#include "../joinedLists/SingleJoinedList.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"

namespace DomainBoundary {
    void apply_boundary_condition(Field &field, SingleJoinedList** index_fields,
                                  BoundaryData* boundary_data, bool sync = true);
}  // namespace DomainBoundary
#endif /* ARTSS_BOUNDARYCONDITION_DOMAINBOUNDARY_H */
