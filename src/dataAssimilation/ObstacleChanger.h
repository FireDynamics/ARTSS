/// \file      ObstacleChanger.h
/// \brief
/// \date      Nov 23, 2022
/// \author    My Linh Wuerzburger
/// \copyright <2015-2022> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_DATAASSIMILATION_OBSTACLECHANGER_H
#define ARTSS_DATAASSIMILATION_OBSTACLECHANGER_H


#include "../interfaces/IParameterReader.h"
#include "../utility/settings/Settings.h"
#include "../utility/Utility.h"
#include "../solver/SolverController.h"

class ObstacleChanger : public IParameterReader {
public:
    ObstacleChanger(const SolverController &solver_controller,
                    const Settings::obstacles_parameters &obstacle_parameters) :
            m_solver_controller(solver_controller),
            m_obstacles(obstacle_parameters),
            m_logger(Utility::create_logger(typeid(this).name())) { }

    return_parameter_reader read_config(const std::string &filename) override;

private:
    const SolverController &m_solver_controller;
    const Settings::obstacles_parameters &m_obstacles;
    std::shared_ptr<spdlog::logger> m_logger;
};


#endif /* ARTSS_DATAASSIMILATION_OBSTACLECHANGER_H */
