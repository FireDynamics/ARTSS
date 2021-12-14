/// \file       BoundaryDataController.h
/// \brief      Controller class for boundary data
/// \date       Dec 09, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_BOUNDARYDATACONTROLLER_H_
#define ARTSS_BOUNDARY_BOUNDARYDATACONTROLLER_H_

#include <vector>

#include "BoundaryData.h"
#include "../field/Field.h"
#include "../joinedLists/SingleJoinedList.h"
#include "../joinedLists/MultipleJoinedList.h"
#include "../utility/Utility.h"
#include "../utility/settings/Settings.h"


class BoundaryDataController {
 public:
    explicit BoundaryDataController(const std::vector<Settings::BoundarySetting> &boundary);
    ~BoundaryDataController() = default;
    void apply_boundary_condition(
            Field &field,
            SingleJoinedList **index_fields,
            bool sync = false);
    void apply_boundary_condition_obstacle(
            Field &field,
            MultipleJoinedList **index_fields,
            size_t id,
            bool sync = false);
    void apply_boundary_condition_surface(
            Field &field,
            MultipleJoinedList **index_fields,
            size_t id,
            bool sync = false);
    void print() const;
    std::vector<FieldType> get_used_fields();

 private:
    std::vector<BoundaryData> m_boundary_data;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    void add_boundary_data(const Settings::BoundarySetting& boundary);
};

#endif /* ARTSS_BOUNDARY_BOUNDARYDATACONTROLLER_H_ */

