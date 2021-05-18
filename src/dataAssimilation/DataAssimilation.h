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
#include "../visualisation/FieldIO.h"

struct AssimilationMethods {
    inline static const std::string Test = "Test";
};

class DataAssimilation {
 public:
    DataAssimilation(FieldController *field_controller);
    void assimilate(real t_cur);
 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    FieldController *m_field_controller;
    bool m_assimilated = false;

    IDataAssimilationFunction *func;
    FieldIO *m_field_IO;
};


#endif /* ARTSS_SRC_DATAASSIMILATION_DATAASSIMILATION_H */
