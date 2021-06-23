/// \file      DataAssimilation.h
/// \brief
/// \date      May 18, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_SRC_DATAASSIMILATION_DATAASSIMILATION_H
#define ARTSS_SRC_DATAASSIMILATION_DATAASSIMILATION_H


#include <string>
#include "../field/FieldController.h"
#include "../interfaces/IDataAssimilationFunction.h"
#include "../utility/Utility.h"
#include "FieldIO.h"

struct AssimilationMethods {
    inline static const std::string Standard = "default";
};

class DataAssimilation {
 public:
    explicit DataAssimilation(const FieldController &field_controller);
    void save_data(real t_cur);

    bool requires_rollback() const { return m_rollback; }

    real get_new_time_value() const;

    void initiate_rollback();

private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    FieldController m_field_controller;
    bool m_assimilated = false;

    IDataAssimilationFunction *m_func;

    bool m_rollback = false;
    real m_t_cur = -1;

    void assimilate(std::string &file_name);
#ifdef ASSIMILATION
    void config_MPI();
#endif
};


#endif /* ARTSS_SRC_DATAASSIMILATION_DATAASSIMILATION_H */
