/// \file       Settings.h
/// \brief      Access parameters of XML file
/// \date       Nov 17, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_SETTINGS_H
#define ARTSS_UTILITY_SETTINGS_H

#include "tinyxml2.h"
#include "BoundarySetting.h"
#include "ObstacleSetting.h"
#include "SurfaceSetting.h"
#include "../GlobalMacrosTypes.h"

#ifndef BENCHMARKING
#include <spdlog/logger.h>
#include <memory>
#endif

#include <map>
#include <string>
#include <vector>
#include <unordered_map>


namespace Settings {
class Settings {
 public:
     explicit Settings(const std::string& path);
     void print_config() const;

     std::string get(std::string path) const;
     std::string sget(std::string path) const;
     void set(std::string path, std::string val) {
         sset(path, val);
#ifndef BENCHMARKING
         m_logger->debug(R"(set: "{}" to value: "{}")", path, val);
#endif
     }
     void sset(std::string path, std::string val) { m_proxy.insert({path, val}); }

     bool get_bool(std::string path) const { return get(path) == "Yes"; }
     int get_int(std::string path) const { return std::stoi(get(path)); }
     int get_size_t(std::string path) const { return std::stol(get(path)); }
     real get_real(std::string path) const { return real(std::stod(get(path))); }

     std::string get_filename() const { return filename; }

     std::vector<BoundarySetting> get_boundaries() const { return m_boundaries; }
     std::vector<ObstacleSetting> get_obstacles() const { return m_obstacles; }
     std::vector<SurfaceSetting> get_surfaces() const { return m_surfaces; }

 private:
     void read_config(std::string prefix, tinyxml2::XMLElement *elem);
     void read_boundaries(tinyxml2::XMLElement *elem);
     void read_obstacles(tinyxml2::XMLElement *elem);
     void read_surfaces(tinyxml2::XMLElement *elem);

     std::vector<BoundarySetting> m_boundaries;
     std::vector<ObstacleSetting> m_obstacles;
     std::vector<SurfaceSetting> m_surfaces;

     std::unordered_multimap<std::string, std::string> m_proxy;
     std::string filename;

#ifndef BENCHMARKING
     std::shared_ptr<spdlog::logger> m_logger;
#endif
};
}
#endif
