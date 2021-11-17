/// \file       Settings.h
/// \brief      Access parameters of XML file
/// \date       Nov 17, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_SETTINGS_H
#define ARTSS_UTILITY_SETTINGS_H

#include <string>
#include <unordered_map>

#include "Parameters.h"
#include "GlobalMacrosTypes.h"

#ifndef BENCHMARKING
#include <memory>
#include "Utility.h"
#endif

union settings_obj {
    bool bool_val;
    int int_val;
    real real_val;
    std::string *str_val;
};

class Settings {
 public:
     Settings() : m_params(Parameters::getInstance())
#ifndef BENCHMARKING
    , m_logger(Utility::create_logger(typeid(this).name()))
#endif
    {}

     std::string get(std::string path);
     void set(std::string path, std::string val) {
         m_proxy.insert({path, val});
#ifndef BENCHMARKING
         m_logger->debug("set: \"{}\" to value: \"{}\"", path, val);
#endif
     }

     bool get_bool(std::string path) {
         return get(path) == "Yes";
     }

     int get_int(std::string path) {
         return real(std::stoi(get(path)));
     }

     real get_real(std::string path) {
         return real(std::stod(get(path)));
     }

 private:
     Parameters *m_params;
     std::unordered_multimap<std::string, std::string> m_proxy;

#ifndef BENCHMARKING
     std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif
