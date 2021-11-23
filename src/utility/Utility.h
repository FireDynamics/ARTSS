/// \file       Utility.h
/// \brief      Offers some tools
/// \date       October 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_UTILITY_UTILITY_H_
#define ARTSS_UTILITY_UTILITY_H_

#include <string>
#include <vector>
#include "Settings.h"
#include "GlobalMacrosTypes.h"

#ifndef BENCHMARKING
#define FMT_USE_UDL_TEMPLATE 0

#include "spdlog/logger.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#endif

class Field;

namespace Utility {
    size_t get_index(real physical_coordinate, real spacing, real start_coordinate);
    std::vector<std::string> split(const char *text, char delimiter);
    std::vector<std::string> split(const std::string &text, char delimiter);
    std::vector<size_t> mergeSortedListsToUniqueList(size_t *list1, size_t size_list1, size_t *list2, size_t size_list2);
    std::string remove_extension(const std::string &filename);
    void log_field_info(Settings const &settings, Field &field, const std::string &text, const std::string &logger_name);

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> create_logger(Settings const &settings, std::string loggerName);
#endif
}  // namespace Utility

#endif /* ARTSS_UTILITY_UTILITY_H_ */

