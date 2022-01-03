/// \file       Utility.h
/// \brief      Offers some tools
/// \date       October 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_UTILITY_UTILITY_H_
#define ARTSS_UTILITY_UTILITY_H_

#include <string>
#include <vector>
#include <tuple>

#include "GlobalMacrosTypes.h"
#include "settings/Settings.h"
#include "../field/Field.h"

#ifndef BENCHMARKING
#define FMT_USE_UDL_TEMPLATE 0

#include "spdlog/logger.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <memory>
#endif


class Field;

namespace Utility {
    std::tuple<size_t, size_t, size_t> get_coordinates(size_t index, size_t Nx, size_t Ny);
    size_t get_index(real physical_coordinate, real spacing, real start_coordinate);
    std::string to_upper(std::string string);
    std::string to_lower(std::string string);
    std::vector<std::string> split(const char *text, char delimiter);
    std::vector<std::string> split(const std::string &text, char delimiter);
    std::string remove_extension(const std::string &filename);
    void log_field_info(Field &field, const std::string &text, const std::string &logger_name);
    std::string trim(const std::string &string);
    void calculate_whole_and_decimal_digits(std::string &number1, std::string &number2, int *decimal_number_digits, int *whole_number_digits);
#ifndef BENCHMARKING
    void create_logger(const std::string &log_level, const std::optional<std::string> &log_file);
    std::shared_ptr<spdlog::logger> create_logger(std::string const& logger_name);
#endif
}  // namespace Utility

#endif /* ARTSS_UTILITY_UTILITY_H_ */

