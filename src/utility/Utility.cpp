/// \file       Utility.cpp
/// \brief      Offers some tools
/// \date       October 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>
#include <sstream>
#include "Utility.h"
#include "GlobalMacrosTypes.h"
#include "Parameters.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"
#include "../field/Field.h"

#ifndef BENCHMARKING
#include <spdlog/cfg/helpers.h>
#include <clocale>
#endif


namespace Utility {
//======================================== get index =====================================
// ***************************************************************************************
/// \brief  Snaps value to grid discretisation
/// \param  physical_coordinate physical coordinate
/// \param  spacing dx/dy/dz
/// \param  start_coordinate X1/Y1/Z1
/// \return real Calculated index (i/j/k) in (x/y/z)-direction
// ***************************************************************************************
size_t get_index(real physical_coordinate, real spacing, real start_coordinate) {
    return std::round((-start_coordinate + physical_coordinate) / spacing) + 1;
}

// ================================= Split string at character ==========================================
// ***************************************************************************************
/// \brief  Splits a string at a defined char
/// \param  text     String
/// \param  delimiter Where to split
// ***************************************************************************************
std::vector<std::string> split(const std::string &text, char delimiter) {
    return split(text.c_str(), delimiter);
}

// ================================= Split string at character ==========================================
// ***************************************************************************************
/// \brief  Splits a string at a defined char
/// \param  text     String
/// \param  delimiter Where to split
// ***************************************************************************************
std::vector<std::string> split(const char *text, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(text);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

#ifndef BENCHMARKING
// ======================= creates a new logger ================================
// *****************************************************************************
/// \brief  creates a new named logger this function is only available if BENCHMARKING is not enabled
/// \param  loggerName name of logger, written to log file
// *****************************************************************************
std::shared_ptr<spdlog::logger> create_logger(std::string logger_name, std::string log_level, std::string log_file) {
    static std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> stdout_sink;
    static std::shared_ptr<spdlog::sinks::basic_file_sink_mt> file_sink;

    if (!stdout_sink) {
        stdout_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        auto level = spdlog::level::from_str(log_level);
        stdout_sink->set_level(level);
        stdout_sink->set_pattern("%^%-8l: %v%$");
    }

    if (!file_sink) {
        file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, false);
        file_sink->set_level(spdlog::level::trace);
    }

    std::vector<spdlog::sink_ptr> sinks;
    sinks.reserve(2);
    sinks.push_back(stdout_sink);
    sinks.push_back(file_sink);
    auto logger = std::make_shared<spdlog::logger>(logger_name, begin(sinks), end(sinks));
    logger->flush_on(spdlog::level::err);
    logger->set_level(spdlog::level::trace);

    return logger;
}

std::shared_ptr<spdlog::logger> create_logger(std::string logger_name) {
    auto params = Parameters::getInstance();
    std::string log_level = params->get("logging/level");
    std::string log_file = params->get("logging/file");
    return create_logger(logger_name, log_level, log_file);
}
#endif

std::vector<size_t> mergeSortedListsToUniqueList(size_t *list1, size_t size_list1, size_t *list2, size_t size_list2) {
    std::vector<size_t> result;
    size_t counter1 = 0;
    size_t counter2 = 0;

    if (list1[counter1] < list2[counter2]) {
        result.push_back(list1[counter1]);
        counter1++;
    } else {
        result.push_back(list2[counter2]);
        counter2++;
    }
    while (counter1 < size_list1 && counter2 < size_list2) {
        if (list1[counter1] == result[result.size()-1]) {
            counter1++;
            continue;
        }
        if (list2[counter2] == result[result.size()-1]) {
            counter2++;
            continue;
        }
        if (list1[counter1] < list2[counter2]) {
            result.push_back(list1[counter1]);
            counter1++;
        } else {
            result.push_back(list2[counter2]);
            counter2++;
        }
    }
    if (counter1 < size_list1) {
        if (list1[counter1] == result[result.size()-1]) {
            counter1++;
        }
        for (size_t c = counter1; c < size_list1; c++) {
            result.push_back(list1[c]);
        }
    } else {
        if (list2[counter2] == result[result.size()-1]) {
            counter2++;
        }
        for (size_t c = counter2; c < size_list2; c++) {
            result.push_back(list2[c]);
        }
    }
    return result;
}

void log_field_info(Field *field, const std::string& text, const std::string& logger_name) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger(logger_name);
#endif
    auto data = field->data;
    auto boundary = BoundaryController::getInstance();
    size_t *inner_list = boundary->get_inner_list_level_joined();
    size_t size_inner_list = boundary->get_size_inner_list();

    size_t idx = inner_list[0];
    real minimum_inner = data[idx];
    real maximum_inner = data[idx];
    real average_inner = data[idx];
    for (size_t i = 1; i < size_inner_list; i++) {
        idx = inner_list[i];
        real value = data[idx];
        if (value < minimum_inner) {
            minimum_inner = value;
        }
        if (value > maximum_inner) {
            maximum_inner = value;
        }
        average_inner += value;
    }
    average_inner /= size_inner_list;
#ifndef BENCHMARKING
    logger->info("minimum inner {}: {}", text, minimum_inner);
    logger->info("maximum inner {}: {}", text, maximum_inner);
    logger->info("average inner {}: {}", text, average_inner);
#endif

    size_t *boundary_list = boundary->get_boundary_list_level_joined();
    size_t size_boundary_list = boundary->get_size_boundary_list();

    idx = boundary_list[0];
    real minimum_boundary = data[idx];
    real maximum_boundary = data[idx];
    real average_boundary = data[idx];
    for (size_t i = 0; i < size_boundary_list; i++) {
        size_t idx = boundary_list[i];
        real value = data[idx];
        if (value < minimum_boundary) {
            minimum_boundary = value;
        }
        if (value > maximum_boundary) {
            maximum_boundary = value;
        }
        average_boundary += value;
    }
    average_boundary /= size_boundary_list;
#ifndef BENCHMARKING
    logger->info("minimum boundary {}: {}", text, minimum_boundary);
    logger->info("maximum boundary {}: {}", text, maximum_boundary);
    logger->info("average boundary {}: {}", text, average_boundary);
#endif

    size_t *obstacle_list = boundary->get_obstacle_list();
    size_t size_obstacle_list = boundary->get_size_obstacle_list();

    if (size_obstacle_list > 0) {
        idx = obstacle_list[0];
        real minimum_obstacle = data[idx];
        real maximum_obstacle = data[idx];
        real average_obstacle = data[idx];
        for (size_t i = 0; i < size_obstacle_list; i++) {
            size_t idx = obstacle_list[i];
            real value = data[idx];
            if (value < minimum_obstacle) {
                minimum_obstacle = value;
            }
            if (value > maximum_obstacle) {
                maximum_obstacle = value;
            }
            average_obstacle += value;
        }
        average_obstacle /= size_obstacle_list;
#ifndef BENCHMARKING
        logger->info("minimum obstacle {}: {}", text, minimum_obstacle);
        logger->info("maximum obstacle {}: {}", text, maximum_obstacle);
        logger->info("average obstacle {}: {}", text, average_obstacle);
#endif
    }
}

}  // namespace Utility
