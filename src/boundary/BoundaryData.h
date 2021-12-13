/// \file       BoundaryData.h
/// \brief      Data class for boundary data
/// \date       Oct 08, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_BOUNDARYDATA_H_
#define ARTSS_BOUNDARY_BOUNDARYDATA_H_

#include <string>
#include <vector>

#include "../boundary/PatchObject.h"
#include "../field/Field.h"
#include "../utility/Utility.h"
#include "../utility/settings/Settings.h"

const size_t number_of_boundary_conditions = 3;
enum BoundaryCondition : int {
    UNKNOWN_CONDITION = -1,
    NEUMANN = 0,
    DIRICHLET = 1,
    PERIODIC = 2
};

class BoundaryData {
 public:
    BoundaryData();
    ~BoundaryData();
    void print();

    static std::string get_boundary_condition_name(BoundaryCondition bc);

    static BoundaryCondition match_boundary_condition(const std::string &string);

    void add_boundary_condition(Patch const &patches,
                                real value,
                                BoundaryCondition const &boundary_condition);
    BoundaryCondition get_boundary_condition(Patch p) const { return m_boundary_conditions[p]; }
    real get_value(Patch p) const { return m_values[p];}

    bool is_empty() const { return !m_has_values; }

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    BoundaryCondition* m_boundary_conditions;
    real* m_values;

    bool m_has_values = false;
};
#endif /* ARTSS_BOUNDARY_BOUNDARYDATA_H_ */

