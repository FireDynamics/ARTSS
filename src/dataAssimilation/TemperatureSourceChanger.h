//
// Created by linh on 05.01.22.
//

#ifndef ARTSS_DATAASSIMILATION_TEMPERATURESOURCECHANGER_H
#define ARTSS_DATAASSIMILATION_TEMPERATURESOURCECHANGER_H

#include "../interfaces/IParameterReader.h"
#include "../utility/settings/Settings.h"

class TemperatureSourceChanger : public IParameterReader {
public:
    Settings::data_assimilation::changes read_config(const std::string &filename);
};


#endif /* ARTSS_DATAASSIMILATION_TEMPERATURESOURCECHANGER_H */
