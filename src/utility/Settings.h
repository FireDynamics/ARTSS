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


class BoundarySetting {
 public:
    explicit BoundarySetting(tinyxml2::XMLElement *xml_element);

    std::string get_field();
    std::string get_patch();
    std::string get_type();
    real get_value();

 private:
    std::string field;
    std::string patch;
    std::string type;
    real value;
};

class ObstacleSetting {
 public:
    explicit ObstacleSetting(tinyxml2::XMLElement *xml_element);

    std::string get_name();
    real get_ox1();
    real get_ox2();
    real get_oy1();
    real get_oy2();
    real get_oz1();
    real get_oz2();

    std::vector<BoundarySetting> get_boundaries();

 private:
    std::string name;
    real ox1, ox2;
    real oy1, oy2;
    real oz1, oz2;
    std::vector<BoundarySetting> m_boundaries;
};

class SurfaceSetting {
 public:
    explicit SurfaceSetting(tinyxml2::XMLElement *xml_element);

    real get_id() const { return id; }
    real get_sx1() const { return sx1; }
    real get_sx2() const { return sx2; }
    real get_sy1() const { return sy1; }
    real get_sy2() const { return sy2; }
    real get_sz1() const { return sz1; }
    real get_sz2() const { return sz2; }

 private:
    int id;
    real sx1, sx2;
    real sy1, sy2;
    real sz1, sz2;
};


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

     std::string get_filename() const {
         return filename;
     }

     std::vector<BoundarySetting> get_boundaries() const;
     std::vector<ObstacleSetting> get_obstacles() const;
     std::vector<SurfaceSetting> get_surfaces() const;

 private:
     void read_config(std::string prefix, tinyxml2::XMLElement *elem);

     std::unordered_multimap<std::string, std::string> m_proxy;
     std::string filename;

#ifndef BENCHMARKING
     std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif
