/// \file       DomainBoundary.h
/// \brief      Applies boundary condition for domain boundary
/// \date       Feb 03, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARYCONDITION_DOMAINBOUNDARY_H
#define ARTSS_BOUNDARYCONDITION_DOMAINBOUNDARY_H

#include <cstdlib>
#include "../utility/GlobalMacrosTypes.h"
#include "../boundary/BoundaryData.h"

class DomainBoundary {
  public:
    void static apply_boundary_condition(real* data_field, size_t** index_fields, const size_t* patch_starts, const size_t* patch_ends, size_t level, BoundaryData* boundary_data, bool sync = true);
  private:
    void static apply_boundary_condition(real* data_field, const size_t* d_patch, size_t patch_start, size_t patch_end, size_t level, int8_t sign_reference_index, size_t reference_index, real value, int8_t sign);
    void static apply_neumann(real* data_field, size_t* d_patch, Patch patch, size_t patch_starts, size_t patch_end, size_t level, real value);
    void static apply_dirichlet(real* data_field, size_t* d_patch, Patch patch, size_t patch_start, size_t patch_end, size_t level, real value);
    void static apply_periodic(real* data_field, size_t* d_patch, Patch patch, size_t patch_start, size_t patch_end, size_t level);
};

#endif /* ARTSS_BOUNDARYCONDITION_DOMAINBOUNDARY_H */
