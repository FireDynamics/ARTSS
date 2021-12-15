/// \file       BoundaryDataController.h
/// \brief      Controll class for boundary data
/// \date       Dec 09, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_BOUNDARYDATACONTROLLER_H_
#define ARTSS_BOUNDARY_BOUNDARYDATACONTROLLER_H_

#include "../utility/Utility.h"
#include "../utility/settings/Settings.h"
#include "BoundaryData.h"
#include "../field/Field.h"
#include "../joinedLists/SingleJoinedList.h"
#include "../joinedLists/MultipleJoinedList.h"

#include <vector>


class BoundaryDataController {
 public:
    explicit BoundaryDataController(Settings::Settings const &settings);
    ~BoundaryDataController();
    void add_boundary_data(Settings::BoundarySetting boundary);
    void apply_boundary_condition(
            Field &field,
            SingleJoinedList **index_fields,
            bool sync = false);
    void apply_boundary_condition_obstacle(
            Field &field,
            MultipleJoinedList **index_fields,
            size_t id,
            bool sync = false);
    void print();
    std::vector<FieldType> get_used_fields();

 private:
    Settings::Settings const &m_settings;
    BoundaryData** m_boundary_data;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_BOUNDARY_BOUNDARYDATACONTROLLER_H_ */

