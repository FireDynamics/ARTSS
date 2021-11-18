/// \file       Settings.h
/// \brief      Access parameters of XML file
/// \date       Nov 17, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_SETTINGS_H
#define ARTSS_UTILITY_SETTINGS_H

#include <map>
#include <string>
#include <iostream>
#include <unordered_map>

#include "tinyxml2.h"
#include "GlobalMacrosTypes.h"

#ifndef BENCHMARKING
#include <memory>
#include "spdlog/logger.h"
#endif


class Settings {
 public:
     explicit Settings(std::string path);
     void print_config() const;

     std::string get(std::string path) const;
     std::string sget(std::string path) const;
     void set(std::string path, std::string val) {
         sset(path, val);
#ifndef BENCHMARKING
         m_logger->debug("set: \"{}\" to value: \"{}\"", path, val);
#endif
     }
     void sset(std::string path, std::string val) {
         m_proxy.insert({path, val});
     }

     bool get_bool(std::string path) const {
         return get(path) == "Yes";
     }

     int get_int(std::string path) const {
         return real(std::stoi(get(path)));
     }

     real get_real(std::string path) const {
         return real(std::stod(get(path)));
     }

 private:
     void read_config(std::string prefix, tinyxml2::XMLElement *elem);

     std::unordered_multimap<std::string, std::string> m_proxy;

#ifndef BENCHMARKING
     std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif
