/// \file      IDataAssimilationWriter.h
/// \brief
/// \date      Aug 25, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#ifndef ARTSS_SRC_INTERFACES_IDATAASSIMILATIONWRITER_H
#define ARTSS_SRC_INTERFACES_IDATAASSIMILATIONWRITER_H

class IDataAssimilationWriter {
public:
    virtual ~IDataAssimilationWriter() = default;
    virtual std::string add_header_data() = 0;
    virtual std::string update_header_data() = 0;
    virtual void parse_header_data(std::string header) = 0;
    virtual std::string add_body_data() = 0;
    virtual std::string update_body_data() = 0;
    virtual void parse_body_data(std::string header) = 0;
};

#endif /* ARTSS_SRC_INTERFACES_IDATAASSIMILATIONWRITER_H */
