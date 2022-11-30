/// \file      ParameterReader.h
/// \brief
/// \date      Jan 03, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_DATAASSIMILATION_PARAMETERREADER_H
#define ARTSS_DATAASSIMILATION_PARAMETERREADER_H

#include <memory>

#include "../interfaces/IParameterReader.h"
#include "../utility/settings/Settings.h"
#include "../utility/Utility.h"
#include "../solver/SolverController.h"

class ParameterReader : public IParameterReader {
public:
    explicit ParameterReader(const SolverController &solver_controller) : m_solver_controller(solver_controller), m_logger(Utility::create_logger(typeid(this).name())) { }

    ~ParameterReader() = default;

    return_parameter_reader read_config(const std::string &file_name) override;

private:
    const SolverController &m_solver_controller;
    std::shared_ptr<spdlog::logger> m_logger;

    bool temperature_source_changer(const tinyxml2::XMLElement *doc, const std::string &context) const;

    bool obstacle_changer(const tinyxml2::XMLElement *head, const std::string &context) const;
};

#endif /* ARTSS_DATAASSIMILATION_PARAMETERREADER_H */
