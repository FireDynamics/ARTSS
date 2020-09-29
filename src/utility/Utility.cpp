/// \file       Utility.cpp
/// \brief      Offers some tools
/// \date       October 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <sstream>
#include "Utility.h"
#include "GlobalMacrosTypes.h"
#include "Parameters.h"

#ifndef BENCHMARKING
#include <spdlog/cfg/helpers.h>
#include <clocale>
#endif


//======================================== get index =====================================
// ***************************************************************************************
/// \brief  Snaps value to grid discretisation
/// \param  physical_coordinate physical coordinate
/// \param  spacing dx/dy/dz
/// \param  start_coordinate X1/Y1/Z1
/// \return real Calculated index (i/j/k) in (x/y/z)-direction
// ***************************************************************************************
size_t Utility::get_index(real physical_coordinate, real spacing, real start_coordinate) {
    return std::round((-start_coordinate + physical_coordinate) / spacing) + 1;
}

// ================================= Calc i,j,k ==========================================
// ***************************************************************************************
/// \brief  calculates indices <i,j,k> from linear index idx
/// \param  idx     linear (global) index
/// \param  Nx      number of cells in x-direction of physical domain
/// \param  Ny      number of cells in y-direction of physical domain
// ***************************************************************************************
std::vector<size_t> Utility::coordinateFromLinearIndex(size_t idx, size_t Nx, size_t Ny) {
    std::vector<size_t> coord;
    size_t k = getCoordinateK(idx, Nx, Ny);
    size_t j = getCoordinateJ(idx, Nx, Ny, k);
    size_t i = getCoordinateI(idx, Nx, Ny, j, k);

    coord.push_back(i);
    coord.push_back(j);
    coord.push_back(k);

    return coord;
}

// ================================= Split string at character ==========================================
// ***************************************************************************************
/// \brief  Splits a string at a defined char
/// \param  text     String
/// \param  delimiter Where to split
// ***************************************************************************************
std::vector<std::string> Utility::split(const std::string &text, char delimiter) {
    return split(text.c_str(), delimiter);
}

// ================================= Split string at character ==========================================
// ***************************************************************************************
/// \brief  Splits a string at a defined char
/// \param  text     String
/// \param  delimiter Where to split
// ***************************************************************************************
std::vector<std::string> Utility::split(const char *text, char delimiter) {
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
/// \brief  creates a new named logger this function is only available if
//          BENCHMARKING is not enabled
/// \param  loggerName name of logger, written to log file
// *****************************************************************************
std::shared_ptr<spdlog::logger> Utility::create_logger(std::string logger_name) {
    static std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> stdout_sink;
    static std::shared_ptr<spdlog::sinks::basic_file_sink_mt> file_sink;

    auto params = Parameters::getInstance();
    std::string log_level = params->get("logging/level");
    std::string log_file = params->get("logging/file");

    if (!stdout_sink) {
        stdout_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        auto level = spdlog::level::from_str(log_level);
        stdout_sink->set_level(level);
        stdout_sink->set_pattern("%^%-8l: %v%$");
    }

    if (!file_sink) {
        file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
            log_file,
            false);
        file_sink->set_level(spdlog::level::trace);
    }

    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(stdout_sink);
    sinks.push_back(file_sink);
    auto logger = std::make_shared<spdlog::logger>(logger_name,
                                                   begin(sinks),
                                                   end(sinks));
    logger->flush_on(spdlog::level::err);
    logger->set_level(spdlog::level::trace);

    return logger;
}
#endif
