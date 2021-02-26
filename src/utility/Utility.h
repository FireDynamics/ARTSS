/// \file       Utility.h
/// \brief      Offers some tools
/// \date       October 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_UTILITY_UTILITY_H_
#define ARTSS_UTILITY_UTILITY_H_

#include <memory>
#include <string>
#include <vector>
#include "GlobalMacrosTypes.h"

#ifndef BENCHMARKING
#define FMT_USE_UDL_TEMPLATE 0

#include "spdlog/logger.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

#endif


namespace Utility {
    size_t get_index(real physical_coordinate, real spacing, real start_coordinate);
    std::vector<std::string> split(const char* text, char delimiter);
    std::vector<std::string> split(const std::string& text, char delimiter);

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> create_logger(std::string loggerName);
    std::shared_ptr<spdlog::logger> create_logger(std::string loggerName, std::string log_level, std::string log_file);
#endif
} // namespace Utility


#endif /* ARTSS_UTILITY_UTILITY_H_ */
