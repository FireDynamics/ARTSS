/// \file       BoundaryData.h
/// \brief      Data class for boundary data
/// \date       Oct 08, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DOMAIN_BOUNDARYDATA_H_
#define ARTSS_DOMAIN_BOUNDARYDATA_H_

#include <string>
#include <optional>
#include <vector>

#include "../domain/PatchObject.h"
#include "../field/Field.h"
#include "../utility/Utility.h"
#include "../utility/settings/Settings.h"

class BoundaryData {
 public:
    explicit BoundaryData();
    ~BoundaryData() = default;
    void print() const;

    void add_boundary_condition(Patch const &patches,
                                std::optional<real> value,
                                BoundaryCondition const &boundary_condition);
    BoundaryCondition get_boundary_condition(Patch p) const { return m_boundary_conditions[p]; }
    real get_value(Patch p) const { return m_values[p].value(); }

    bool is_empty() const { return !m_has_values; }

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    std::vector<BoundaryCondition> m_boundary_conditions;
    std::vector<std::optional<real>> m_values;

    bool m_has_values = false;
};
#endif /* ARTSS_DOMAIN_BOUNDARYDATA_H_ */

