/// \file      DataAssimilation.h
/// \brief
/// \date      May 18, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_SRC_DATAASSIMILATION_DATAASSIMILATION_H
#define ARTSS_SRC_DATAASSIMILATION_DATAASSIMILATION_H


#include <string>
#include "../field/FieldController.h"
#include "../utility/Utility.h"
#include "FieldIO.h"
#include "../solver/SolverController.h"
#include "ParameterReader.h"

struct AssimilationMethods {
    inline static const std::string None = "default";
    inline static const std::string HRRChanger = "HRRChanger";
};

class DataAssimilation {
 public:
    DataAssimilation(const SolverController &solver_controller,
                     FieldController *field_controller,
                     const Settings::assimilation_parameters &settings);
    void save_data(real t_cur);

    bool requires_rollback();
    void initiate_rollback();

    real get_new_time_value() const;

private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    const Settings::assimilation_parameters &m_settings;
    FieldController *m_field_controller;
    const SolverController &m_solver_controller;

    FieldIO *m_field_IO_handler;
    ParameterReader *m_parameter_handler;

    real m_t_cur = -1;

    void read_new_data(std::string &file_name);
    void config_rollback(const char *msg);

    Field m_new_field_u, m_new_field_v, m_new_field_w, m_new_field_p, \
          m_new_field_T, m_new_field_C;
};


#endif /* ARTSS_SRC_DATAASSIMILATION_DATAASSIMILATION_H */
