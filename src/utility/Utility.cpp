/// \file       Utility.cpp
/// \brief      Offers some tools
/// \date       October 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>
#include <sstream>
#include "Utility.h"
#include "GlobalMacrosTypes.h"
#include "../DomainData.h"
#include "../boundary/DomainController.h"

#ifndef BENCHMARKING
#include <spdlog/spdlog.h>
#include <spdlog/cfg/helpers.h>
#include <clocale>
#endif


namespace Utility {
    const static std::string class_name = "Utility";
    const static std::string global_logger = "ARTSS";

    // do not use only for debug purpose
std::vector<size_t> get_coordinates(size_t index, size_t Nx, size_t Ny) {
    std::vector<size_t> coordinates;
    coordinates.reserve(3);
    size_t k = getCoordinateK(index, Nx, Ny);
    size_t j = getCoordinateJ(index, Nx, Ny, k);
    size_t i = getCoordinateI(index, Nx, Ny, j, k);
    coordinates.push_back(i);
    coordinates.push_back(j);
    coordinates.push_back(k);
    return coordinates;
}

//======================================== get index ===============================================
// *************************************************************************************************
/// \brief  Snaps value to grid discretisation
/// \param  physical_coordinate physical coordinate
/// \param  spacing dx/dy/dz
/// \param  start_coordinate X1/Y1/Z1
/// \return real Calculated index (i/j/k) in (x/y/z)-direction
// *************************************************************************************************
size_t get_index(real physical_coordinate, real spacing, real start_coordinate) {
    return std::round((-start_coordinate + physical_coordinate) / spacing) + 1;
}

// ================================= Split string at character =====================================
// *************************************************************************************************
/// \brief  Splits a string at a defined char
/// \param  text     String
/// \param  delimiter Where to split
// *************************************************************************************************
std::vector<std::string> split(const std::string &text, char delimiter) {
    return split(text.c_str(), delimiter);
}

// ================================= Split string at character =====================================
// *************************************************************************************************
/// \brief  Splits a string at a defined char
/// \param  text     String
/// \param  delimiter Where to split
// *************************************************************************************************
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
std::shared_ptr<spdlog::logger> create_logger(const std::string &logger_name) {
    auto &sinks = spdlog::get(global_logger)->sinks();
    auto logger = std::make_shared<spdlog::logger>(logger_name, sinks.begin(), sinks.end());
    logger->sinks()[0]->set_level(sinks[0]->level());
    logger->sinks()[1]->set_level(sinks[1]->level());
    logger->set_level(sinks[0]->level());
    return logger;
}

// ======================= creates a new logger ================================
// *****************************************************************************
/// \brief  creates a new named logger this function is only available
///         if BENCHMARKING is not enabled
/// \param  settings the settings to create the logger
///         ("logging/level", "logging/file")
/// \param  logger_name name of logger, written to log file
// *****************************************************************************
void create_logger(const std::string &log_level, const std::optional<std::string> &log_file) {
    std::vector<spdlog::sink_ptr> sinks;
    static auto stdout_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    auto level = spdlog::level::from_str(log_level);
    stdout_sink->set_level(level);
    stdout_sink->set_pattern("%^%-8l: %v%$");
    sinks.emplace_back(stdout_sink);

    static auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file.value_or("file.log"), false);
    file_sink->set_level(spdlog::level::trace);
    sinks.emplace_back(file_sink);

    auto combined_logger = std::make_shared<spdlog::logger>(global_logger, begin(sinks), end(sinks));
    spdlog::register_logger(combined_logger);  // needed if spdlog::get() should be used elsewhere
    combined_logger->flush_on(spdlog::level::err);
    combined_logger->set_level(spdlog::level::trace);
}
#endif


void log_field_info(Field &field, const std::string &text, const std::string &logger_name) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger(logger_name);
#endif
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    size_t idx = domain_inner_list[0];
    real minimum_inner = field[idx];
    real maximum_inner = field[idx];
    real average_inner = field[idx];
    for (size_t i = 1; i < size_domain_inner_list; i++) {
        idx = domain_inner_list[i];
        real value = field[idx];
        if (value < minimum_inner) {
            minimum_inner = value;
        }
        if (value > maximum_inner) {
            maximum_inner = value;
        }
        average_inner += value;
    }
    average_inner /= static_cast<real>(size_domain_inner_list);
#ifndef BENCHMARKING
    logger->info("minimum inner {}: {}", text, minimum_inner);
    logger->info("maximum inner {}: {}", text, maximum_inner);
    logger->info("average inner {}: {}", text, average_inner);
#endif
}

//================================= Remove extension ===============================================
// *************************************************************************************************
/// \brief  Removes extension from filename
/// \param  filename    xml-file name (via argument)
// *************************************************************************************************
std::string remove_extension(const std::string &filename) {
    size_t lastdot = filename.find_last_of('.');
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}

}  // namespace Utility
