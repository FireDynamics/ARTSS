/// \file 		ObstacleBoundary.h
/// \brief 		Applies boundary condition for obstacle boundary
/// \date 		Feb 03, 2020
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARYCONDITION_OBSTACLEBOUNDARY_H
#define ARTSS_BOUNDARYCONDITION_OBSTACLEBOUNDARY_H

#include <cstdlib>
#include "../Utility/GlobalMacrosTypes.h"
#include "../Boundary/BoundaryData.h"

class ObstacleBoundary {
public:
    void static applyBoundaryCondition(real* data, size_t** indexFields, const size_t* patch_start, const size_t* patch_end, size_t level, BoundaryData* boundaryData, size_t id, bool sync=true);
private:
    void static applyBoundaryCondition(real* dataField, const size_t* d_patch, size_t patch_start, size_t patch_end, size_t level, int referenceIndex, real value, int8_t sign);
    void static applyNeumann(real* dataField, size_t* d_patch, Patch p, size_t patch_start, size_t patch_end, size_t level, real value);
    void static applyDirichlet(real* dataField, size_t* d_patch, Patch p, size_t patch_start, size_t patch_end, size_t level, real value);
    void static applyPeriodic(real* dataField, size_t* d_patch, Patch p, size_t patch_start, size_t patch_end, size_t level, size_t id);
};

#endif /* ARTSS_BOUNDARYCONDITION_OBSTACLEBOUNDARY_H */
