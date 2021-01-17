/// \file       Parameters.cpp
/// \brief      Access parameters of XML file
/// \date       May 20, 2016
/// \author     Arnold
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#include "Parameters.h"
#include "Utility.h"


Parameters *Parameters::single = nullptr;

// Singleton
Parameters *Parameters::getInstance() {
    if (single == nullptr) {
        single = new Parameters();
    }
    return single;
}

// ==================================== Parse ==================================
// *****************************************************************************
/// \brief  parses xml file
/// \param  filename        string (name of xml-file)
// *****************************************************************************
void Parameters::parse(const std::string& filename) {
#ifndef BENCHMARKING
    auto inpFileLogger = Utility::create_logger("XMLFile");
#endif

    if (filename.empty()) {
#ifndef BENCHMARKING
        inpFileLogger->error("no XML file specified, skip reading parameter");
#endif
        return;
    }

    tinyxml2::XMLError eResult = this->doc->LoadFile(filename.c_str());
    m_filename = filename;

#ifndef BENCHMARKING
    std::ifstream inpFileStream(filename, std::ifstream::in);
    if (inpFileStream.is_open()) {
        std::string line;
        while (std::getline(inpFileStream, line)) {
            inpFileLogger->debug(line);
        }

        inpFileStream.close();
    } else {
        inpFileLogger->error("Can't open file: {}", filename);
    }

    inpFileLogger->debug("<!-- eResult={}, filename={} -->", eResult, filename);
#endif
}

// =================================== Getter ==================================
// *****************************************************************************
/// \brief  gets raw string (from xml-file)
/// \param  raw_path   tree path (as string) of xml-file
// *****************************************************************************
std::string Parameters::get(const std::string &raw_path) {
    auto path = Utility::split(raw_path, '/');

    auto last_path_element = path.back();
    path.pop_back();

    auto current = this->doc->RootElement();

    for (auto const& cpath : path) {
        current = current->FirstChildElement(cpath.c_str());
    }

    if (current->Attribute(last_path_element.c_str())) {
        return std::string(current->Attribute(last_path_element.c_str()));
    } else {
        current = current->FirstChildElement(last_path_element.c_str());
        return std::string(current->GetText());
    }
}

// ***************************************************************************************
/// \brief  gets real number (from xml-file)
/// \param  raw_path        tree path (as string) of xml-file
// ***************************************************************************************
real Parameters::get_real(const std::string &raw_path) {
    auto raw_result = this->get(raw_path);
    return real(std::stod(raw_result));
}

// ***************************************************************************************
/// \brief  gets double number (from xml-file)
/// \param  raw_path        tree path (as string) of xml-file
// ***************************************************************************************
double Parameters::get_double(const std::string &raw_path) {
    auto raw_result = this->get(raw_path);
    return std::stod(raw_result);
}

// ***************************************************************************************
/// \brief  gets integer number (from xml-file)
/// \param  raw_path        tree path (as string) of xml-file
// ***************************************************************************************
int Parameters::get_int(const std::string &raw_path) {
    auto raw_result = this->get(raw_path);
    return std::stoi(raw_result);
}

tinyxml2::XMLElement* Parameters::get_first_child(const char* raw_path){
    return this->doc->RootElement()->FirstChildElement(raw_path);
}
