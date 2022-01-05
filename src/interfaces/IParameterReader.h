/// \file      IParameterReader.h
/// \brief
/// \date      Jan 05, 2022
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_INTERFACES_IPARAMETERREADER_H
#define ARTSS_INTERFACES_IPARAMETERREADER_H

#include <string>

class IParameterReader {
    virtual std::string read_config(const std::string &filename) = 0;
};

#endif /* ARTSS_INTERFACES_IPARAMETERREADER_H */