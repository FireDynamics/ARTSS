/// \file      IParameterReader.h
/// \brief
/// \date      Jan 05, 2022
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_INTERFACES_IPARAMETERREADER_H
#define ARTSS_INTERFACES_IPARAMETERREADER_H

#include <string>
#include "../utility/settings/Settings.h"

class IParameterReader {
public:
    virtual Settings::data_assimilation::field_changes read_config(const std::string &filename, real t_cur) = 0;
};

#endif /* ARTSS_INTERFACES_IPARAMETERREADER_H */
