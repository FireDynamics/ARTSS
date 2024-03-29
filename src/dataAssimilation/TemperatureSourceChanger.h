/// \file       TemperatureSourceChanger.h
/// \brief      Class for changing temperature source
/// \date       Jan 05, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2022> Forschungszentrum Juelich All rights reserved.

#ifndef ARTSS_DATAASSIMILATION_TEMPERATURESOURCECHANGER_H
#define ARTSS_DATAASSIMILATION_TEMPERATURESOURCECHANGER_H

#include <stdexcept>

#include "../solver/SolverController.h"
#include "../interfaces/IParameterReader.h"
#include "../utility/settings/Settings.h"
#include "../utility/Utility.h"

class TemperatureSourceChanger : public IParameterReader {
 public:
    TemperatureSourceChanger(const SolverController &solver_controller,
                             const Settings::solver::temperature_source &temperature_source) :
            m_solver_controller(solver_controller),
            m_temperature_source(temperature_source),
            m_logger(Utility::create_logger(typeid(this).name())) { }
    return_parameter_reader read_config(const std::string &filename) override;
 private:
    const SolverController &m_solver_controller;
    const Settings::solver::temperature_source &m_temperature_source;
    std::shared_ptr<spdlog::logger> m_logger;
};

#endif /* ARTSS_DATAASSIMILATION_TEMPERATURESOURCECHANGER_H */
