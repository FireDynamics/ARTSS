/// \file      HRRChanger.cpp
/// \brief
/// \date      Aug 04, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include "HRRChanger.h"
#include "../utility/Parameters.h"
#include "../Domain.h"

HRRChanger::HRRChanger(ISourceFunction *source_function) : m_source_function(source_function) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
}

std::string HRRChanger::add_header_data() {
    return std::string();
}

std::string HRRChanger::update_header_data() {
    return std::string();
}

void HRRChanger::parse_header_data(std::string header) {

}

std::string HRRChanger::add_body_data() {
    return std::string();
}

std::string HRRChanger::update_body_data() {
    return std::string();
}

void HRRChanger::parse_body_data(std::string header) {

}
