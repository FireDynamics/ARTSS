/// \file       Utility.h
/// \brief      Offers some tools
/// \date       October 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_UTILITY_UTILITY_H_
#define ARTSS_UTILITY_UTILITY_H_

#include <string>
#include <vector>

#include "GlobalMacrosTypes.h"
#include "settings/Settings.h"

#ifndef BENCHMARKING
#ifdef _OPENACC
#define GPU_DEBUG
#endif
#define FMT_USE_UDL_TEMPLATE 0

#include "spdlog/logger.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <memory>
#endif


class Field;

namespace Utility {
    std::vector<size_t> get_coordinates(size_t index, size_t Nx, size_t Ny);
    size_t get_index(real physical_coordinate, real spacing, real start_coordinate);
    std::vector<std::string> split(const char *text, char delimiter);
    std::vector<std::string> split(const std::string &text, char delimiter);
    std::string remove_extension(const std::string &filename);
    void log_field_info(Settings::Settings const &settings, Field &field, const std::string &text, const std::string &logger_name);

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> create_logger(
            std::string const logger_level,
            std::string const logger_file,
            std::string const logger_name);
    std::shared_ptr<spdlog::logger> create_logger(
            Settings::Settings const &settings,
            std::string const logger_name);
    std::shared_ptr<spdlog::logger> create_logger(std::string const logger_name);
#endif
#ifdef GPU_DEBUG
    std::shared_ptr<spdlog::logger> create_gpu_logger(std::string logger_name);
    std::shared_ptr<spdlog::logger> create_gpu_logger(Settings::Settings const &settings, std::string logger_name);
#endif
}  // namespace Utility

#endif /* ARTSS_UTILITY_UTILITY_H_ */

