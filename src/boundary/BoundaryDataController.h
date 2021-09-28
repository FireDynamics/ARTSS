/// \file       BoundaryDataController.h
/// \brief      Controll class for boundary data
/// \date       Dec 09, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_BOUNDARYDATACONTROLLER_H_
#define ARTSS_BOUNDARY_BOUNDARYDATACONTROLLER_H_

#include <vector>
#include "../utility/tinyxml2.h"
#include "../utility/Utility.h"
#include "BoundaryData.h"

class BoundaryDataController {
 public:
    BoundaryDataController();
    ~BoundaryDataController();
    void add_boundary_data(tinyxml2::XMLElement *xml_element);
    void apply_boundary_condition(
            Field &data,
            size_t **index_fields,
            const size_t *patch_start, const size_t *patch_end,
            bool sync = false);
    void apply_boundary_condition_obstacle(
            real *data,
            size_t **index_fields,
            size_t *patch_start, size_t *patch_end,
            FieldType field_type,
            size_t level,
            size_t id,
            bool sync = false);
    void print();
    std::vector<FieldType> get_used_fields();

 private:
    BoundaryData** m_boundary_data;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_BOUNDARY_BOUNDARYDATACONTROLLER_H_ */

