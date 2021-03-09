/// \file       Parameters.h
/// \brief      Access parameters of XML file
/// \date       May 20, 2016
/// \author     Arnold
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_PARAMETERS_H
#define ARTSS_UTILITY_PARAMETERS_H

#include <string>
#include <iostream>
#include <fstream>

#include "tinyxml2.h"
#include "GlobalMacrosTypes.h"
#include "Utility.h"


class Parameters {
 public:
    static Parameters* getInstance();
    void parse(const std::string& filename);

    // Getter
    std::string get(const std::string& raw_path);
    real get_real(const std::string& raw_path);
    double get_double(const std::string& raw_path);
    int get_int(const std::string& raw_path);
    std::string get_filename() {return m_filename; }

    void printAllXMLAttributes(std::string prefix, tinyxml2::XMLElement *node);

    tinyxml2::XMLElement *get_first_child(const std::string &raw_path);
    tinyxml2::XMLElement *get_first_child(const char *raw_path);

 private:
    tinyxml2::XMLDocument* doc;
    static Parameters* single;

    Parameters() {this->doc = new tinyxml2::XMLDocument;}
    std::string m_filename;

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_UTILITY_PARAMETERS_H */
