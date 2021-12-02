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
#include "../boundary/BoundaryController.h"
#include "../field/Field.h"

#ifndef BENCHMARKING

#include <spdlog/cfg/helpers.h>
#include <clocale>

#endif


namespace Utility {
    static std::string class_name = "Utility";


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

#ifdef GPU_DEBUG
    // =====================creates a new logger for the GPU =======================
    // *****************************************************************************
    /// \brief  creates a new named logger this function is only available
    ///         if BENCHMARKING is not enabled and _OPENACC is enabled
    /// \param  loggerName name of logger, written to log file
    // *****************************************************************************
        std::shared_ptr<spdlog::logger> create_gpu_logger(std::string logger_name) {
        /*
            static std::shared_ptr<spdlog::sinks::basic_file_sink_mt> file_sink;
            std::cout << "print name: " << logger_name << std::endl;
            auto params = Parameters::getInstance();
            std::string log_level = "debug";
            std::string log_file = params->get("logging/file") + "_gpu";

            if (!file_sink) {
                file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, false);
                file_sink->set_level(spdlog::level::trace);
            }

            std::shared_ptr<spdlog::logger> logger = spdlog::basic_logger_mt(logger_name, log_file);
            logger->flush_on(spdlog::level::err);
            logger->set_level(spdlog::level::trace);
            return logger;
        */

        static std::shared_ptr<spdlog::sinks::basic_file_sink_mt> file_sink;

        auto params = Parameters::getInstance();
        std::string log_level = params->get("logging/level");
        std::string log_file = params->get("logging/file") + "_gpu";

        if (!file_sink) {
            file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, false);
            file_sink->set_level(spdlog::level::trace);
        }

        std::vector<spdlog::sink_ptr> sinks;
        sinks.reserve(1);
        sinks.push_back(file_sink);
        auto logger = std::make_shared<spdlog::logger>(logger_name, begin(sinks), end(sinks));
        logger->flush_on(spdlog::level::err);
        logger->set_level(spdlog::level::trace);
        return logger;
        }
#endif

#ifndef BENCHMARKING
// ======================= creates a new logger ================================
// *****************************************************************************
/// \brief  creates a new named logger this function is only available
///         if BENCHMARKING is not enabled
/// \param  settings the settings to create the logger
///         ("logging/level", "logging/file")
/// \param  loggerName name of logger, written to log file
// *****************************************************************************
std::shared_ptr<spdlog::logger> create_logger(Settings::Settings const &settings, std::string const logger_name) {
    return create_logger(
            settings.sget("logging/level"),
            settings.sget("logging/file"),
            logger_name);
}

// ======================= creates a new logger ================================
// *****************************************************************************
/// \brief  creates a new named logger this function is only available
///         if BENCHMARKING is not enabled
/// \param  level the level of visable messages
/// \param  file the file to write into
/// \param  loggerName name of logger, written to log file
// *****************************************************************************
std::shared_ptr<spdlog::logger> create_logger(
        std::string const log_level,
        std::string const log_file,
        std::string const logger_name) {
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
#endif


void log_field_info(Settings::Settings const &settings, Field &field, const std::string &text, const std::string &logger_name) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger(settings, logger_name);
#endif
    auto boundary = BoundaryController::getInstance();
    size_t *inner_list = boundary->get_domain_inner_list_level_joined();
    size_t size_inner_list = boundary->get_size_domain_inner_list_level_joined(0);

    size_t idx = inner_list[0];
    real minimum_inner = field[idx];
    real maximum_inner = field[idx];
    real average_inner = field[idx];
    for (size_t i = 1; i < size_inner_list; i++) {
        idx = inner_list[i];
        real value = field[idx];
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
