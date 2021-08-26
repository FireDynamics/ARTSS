/// \file      Zero.h
/// \brief
/// \date      Aug 25, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_SRC_DATAASSIMILATION_ZERO_H
#define ARTSS_SRC_DATAASSIMILATION_ZERO_H

#include "../interfaces/IDataAssimilationWriter.h"

class Zero : public IDataAssimilationWriter {
 public:
    Zero() = default;
    std::string add_header_data() override { return {}; }
    std::string update_header_data() override { return {}; }
    void parse_header_data(std::string header) override { }
    std::string add_body_data() override { return {}; }
    std::string update_body_data() override { return {}; }
    void parse_body_data(std::string header) override { }
};

#endif /* ARTSS_SRC_DATAASSIMILATION_ZERO_H */
