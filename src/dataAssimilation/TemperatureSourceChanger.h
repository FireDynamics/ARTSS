//
// Created by linh on 05.01.22.
//

#ifndef ARTSS_DATAASSIMILATION_TEMPERATURESOURCECHANGER_H
#define ARTSS_DATAASSIMILATION_TEMPERATURESOURCECHANGER_H

#include "../solver/SolverController.h"
#include "../interfaces/IParameterReader.h"
#include "../utility/settings/Settings.h"

class TemperatureSourceChanger : public IParameterReader {
public:
    TemperatureSourceChanger(const SolverController &solver_controller,
                             const Settings::solver::temperature_source &temperature_source) :
                             m_solver_controller(solver_controller),
                             m_temperature_source(temperature_source) {};
    Settings::data_assimilation::field_changes read_config(const std::string &filename) override;
private:
    const SolverController &m_solver_controller;
    const Settings::solver::temperature_source &m_temperature_source;
};


#endif /* ARTSS_DATAASSIMILATION_TEMPERATURESOURCECHANGER_H */
