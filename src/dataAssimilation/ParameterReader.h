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

class ParameterReader : public IParameterReader {
 public:
    ParameterReader() : m_logger(Utility::create_logger(typeid(this).name())) {}
    ~ParameterReader() = default;

    return_parameter_reader read_config(const std::string &file_name) override;
 private:
    std::shared_ptr<spdlog::logger> m_logger;
};


#endif /* ARTSS_DATAASSIMILATION_PARAMETERREADER_H */
