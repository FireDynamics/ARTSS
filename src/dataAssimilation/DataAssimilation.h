/// \file      DataAssimilation.h
/// \brief
/// \date      May 18, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_DATAASSIMILATION_DATAASSIMILATION_H
#define ARTSS_DATAASSIMILATION_DATAASSIMILATION_H


#include <memory>
#include <string>
#include "../field/FieldController.h"
#include "../utility/Utility.h"
#include "FieldIO.h"
#include "../solver/SolverController.h"
#include "ParameterReader.h"
#include "../interfaces/IParameterReader.h"

struct AssimilationMethods {
    inline static const std::string standard = "default";
    inline static const std::string temperature_source = "TemperatureSourceChanger";
};

struct DataAssimilationPackageHeader {
    double time;
    int file_name_len;
};

class DataAssimilation {
 public:
    DataAssimilation(const SolverController &solver_controller,
                     FieldController *field_controller,
                     const Settings::Settings &settings);
    void save_data(real t_cur);

    bool requires_rollback(real t_cur);
    void initiate_time_skip(real t_cur);

    real get_new_time_value() const;

    bool load_data();

 private:
    std::shared_ptr<spdlog::logger> m_logger;
    const Settings::Settings &m_settings;
    FieldController *m_field_controller;
    const SolverController &m_solver_controller;

    FieldIO *m_field_IO_handler;
    IParameterReader *m_parameter_handler;

    real m_t_cur = 0;
    real m_output_time_interval;
    int m_time_interval_counter;

    bool config_rollback(const char *msg);

    Field m_new_field_u, m_new_field_v, m_new_field_w, \
          m_new_field_p, m_new_field_T, m_new_field_C;
};

#endif /* ARTSS_DATAASSIMILATION_DATAASSIMILATION_H */
