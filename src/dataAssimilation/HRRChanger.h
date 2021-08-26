/// \file      HRRChanger.h
/// \brief
/// \date      Aug 04, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_SRC_DATAASSIMILATION_HRRCHANGER_H
#define ARTSS_SRC_DATAASSIMILATION_HRRCHANGER_H


#include <string>
#include "FieldIO.h"
#include "../utility/Utility.h"
#include "../interfaces/ISourceFunction.h"

class HRRChanger: public IDataAssimilationWriter {
 public:
    explicit HRRChanger(ISourceFunction *source_function);
    std::string add_header_data() override;
    std::string update_header_data() override;
    void parse_header_data(std::string header) override;
    std::string add_body_data() override;
    std::string update_body_data() override;
    void parse_body_data(std::string header) override;
 private:
    ISourceFunction *m_source_function;

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};


#endif /* ARTSS_SRC_DATAASSIMILATION_HRRCHANGER_H */
