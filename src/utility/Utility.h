/// \file       Utility.h
/// \brief      Offers some tools
/// \date       October 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_UTILITY_UTILITY_H_
#define ARTSS_UTILITY_UTILITY_H_

#include <string>
#include <vector>

#ifndef BENCHMARKING
#define FMT_USE_UDL_TEMPLATE 0
#include "spdlog/logger.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#endif

class Utility {
 public:
    static std::vector<size_t> coordinateFromLinearIndex(size_t idx,
                                                         size_t Nx,
                                                         size_t Ny);
//  static size_t getCoordinateI(size_t idx,
//                               size_t Nx,
//                               size_t Ny,
//                               size_t j,
//                               size_t k);
//  static size_t getCoordinateJ(size_t idx, size_t Nx, size_t Ny, size_t k);
//  static size_t getCoordinateK(size_t idx, size_t Nx, size_t Ny);
    static std::vector<std::string> split(const char* text, char delimiter);
    static std::vector<std::string> split(const std::string& text,
                                          char delimiter);

#ifndef BENCHMARKING
    static std::shared_ptr<spdlog::logger> create_logger(std::string loggerName);
#endif

 private:
    Utility() = default;
};


#endif /* ARTSS_UTILITY_UTILITY_H_ */
