/// \file 		Parameters.h
/// \brief 		Access parameters of XML file
/// \date 		May 20, 2016
/// \author 	Arnold
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_PARAMETERS_H
#define ARTSS_UTILITY_PARAMETERS_H

#include <string>

#include "tinyxml2.h"
#include "GlobalMacrosTypes.h"

class Parameters {

private:
    tinyxml2::XMLDocument* doc;
    static Parameters* single;

    Parameters() {this->doc = new tinyxml2::XMLDocument;};

public:
    static Parameters* getInstance();
    void parse(const std::string& filename);

    // Getter
    std::string get(const std::string& raw_path);
    real getReal(const std::string& raw_path);
    double getDouble(const std::string& raw_path);
    int getInt(const std::string& raw_path);

    tinyxml2::XMLElement *getRootElement() {return doc->RootElement();};
};

#endif /* ARTSS_UTILITY_PARAMETERS_H */
