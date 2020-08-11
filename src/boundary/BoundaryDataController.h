/// \file 		BoundaryDataController.h
/// \brief 		Controll class for boundary data
/// \date 		Dec 09, 2019
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDRY_BOUNDARYDATACONTROLLER_H
#define ARTSS_BOUNDRY_BOUNDARYDATACONTROLLER_H


#include "../utility/tinyxml2.h"
#include "BoundaryData.h"

class BoundaryDataController {
public:
    BoundaryDataController();
    ~BoundaryDataController();
    void addBoundaryData(tinyxml2::XMLElement *xmlElement);
    void applyBoundaryCondition(real *data, size_t **indexFields, size_t *patch_start, size_t *patch_end, FieldType fieldType, size_t level, bool sync=false);
    void applyBoundaryConditionObstacle(real *data, size_t **indexFields, size_t *patch_start, size_t *patch_end, FieldType fieldType, size_t level, size_t id, bool sync=false);
    //void applyBoundaryCondition(real *data, FieldType fieldType, size_t level, bool sync=false);
    void print();

    //void setIndexFields(size_t** indexFields);

    std::vector<FieldType> get_used_fields();

private:
    BoundaryData** m_boundaryData;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

};


#endif /* ARTSS_BOUNDRY_BOUNDARYDATACONTROLLER_H */
