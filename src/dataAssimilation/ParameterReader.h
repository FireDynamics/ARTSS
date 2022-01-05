/// \file      ParameterReader.h
/// \brief
/// \date      Jan 03, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_DATAASSIMILATION_PARAMETERREADER_H
#define ARTSS_DATAASSIMILATION_PARAMETERREADER_H

#include "../interfaces/IParameterReader.h"

class ParameterReader : public IParameterReader {
public:
    ParameterReader() = default;
    ~ParameterReader() = default;
};


#endif /* ARTSS_DATAASSIMILATION_PARAMETERREADER_H */
