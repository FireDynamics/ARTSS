/// \file       main.cpp
/// \brief      Hybrid (CPU and GPU) Version for solving Navier-Stokes equations
///             in 3D
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

#include "interfaces/SolverI.h"
#include "solver/DiffusionTurbSolver.h"
#include "solver/DiffusionSolver.h"
#include "solver/AdvectionSolver.h"
#include "solver/AdvectionDiffusionSolver.h"
#include "solver/NSSolver.h"
#include "solver/NSTurbSolver.h"
#include "solver/NSTempSolver.h"
#include "solver/NSTempTurbSolver.h"
#include "solver/NSTempConSolver.h"
#include "solver/NSTempTurbConSolver.h"
#include "solver/PressureSolver.h"
#include "TimeIntegration.h"
#include "utility/tinyxml2.h"
#include "utility/Parameters.h"

#ifndef PROFILING
#include "analysis/Analysis.h"
#include "utility/Visual.h"
#endif

#ifdef _OPENACC
    #include <openacc.h>
#endif


//================================= parse_params =============================
// ***************************************************************************
/// \brief parses arguments provided to the program for logging
/// \param argc         number of arguments
/// \param argv         argument vector with strings
// ***************************************************************************
void parse_params(int argc, char **argv) {
    char parm;
    std::string log_file, log_level, XMLfilename;
    auto params = Parameters::getInstance();

    while ((parm = static_cast<char>(getopt(argc, argv, "o:l:"))) != -1)
        switch (parm) {
            case 'l':  // loglevel
                log_level.assign(optarg);
                break;

            case 'o':  // output for logging
                log_file.assign(optarg);
                break;

            case '?':
                if (optopt == 'l')  // if -l expects a loglevel
                    spdlog::error("Option -{} requires an argument",
                                  static_cast<char>(optopt));

                else if (optopt == 'o')  // if -o expects a path
                    spdlog::error("Option -{} requires an argument",
                                  static_cast<char>(optopt));

                else  // yet unknown option requested
                    spdlog::error("Unknown option {}",
                                  static_cast<char>(optopt));

                break;

            default:  // we should never reach this point
                spdlog::error("Error in argument parsing");
        }

    // get xml path as last free argument
    if (optind == argc-1)
        XMLfilename.assign(argv[optind]);
    else
        spdlog::error("Filepath is missing");

    spdlog::info("Provided XML-path: {}", XMLfilename);
    spdlog::info("Provided logging level: {}", log_level);
    spdlog::info("Provided logging output: {}", log_file);

    // setting up logging level
    if (log_level != "") {
        auto level = spdlog::level::from_str(log_level);
        spdlog::default_logger()->set_level(level);
    }
    auto def_logger = spdlog::default_logger();
    auto log_level_str = spdlog::level::to_string_view(def_logger->level());
    spdlog::info("We are at log-level: {}", log_level_str);

    // setting up loger output
    if (log_file == "-") {  // ouptut on stdout
    } else {  // output in file
        if (log_file == "")  // default file if not otherwise specified
            log_file = "./log.txt";

        auto file_logger = spdlog::basic_logger_mt("basic_logger", log_file);
        spdlog::set_default_logger(file_logger);
        // make sure we see the last message
        spdlog::default_logger()->flush_on(spdlog::level::err);
    }

    spdlog::info("Parsing of {}", XMLfilename);
    params->parse(XMLfilename);
    spdlog::info("Parsing of {} completed", XMLfilename);
}

int main(int argc, char **argv) {
    // 0. Initialization
    // Parameters
    std::string XMLfilename;
    auto params = Parameters::getInstance();

#ifndef PROFILING
    parse_params(argc, argv);
#else
    if (argc > 1) XMLfilename.assign(argv[1]);
    params->parse(XMLfilename);
#endif

    // Solver
    SolverI* solver;
    std::string string_solver = params->get("solver/description");
    if (string_solver == SolverTypes::DiffusionSolver) {
        solver = new DiffusionSolver();
    } else if (string_solver == SolverTypes::AdvectionSolver) {
        solver = new AdvectionSolver();
    } else if (string_solver == SolverTypes::AdvectionDiffusionSolver) {
        solver = new AdvectionDiffusionSolver();
    } else if (string_solver == SolverTypes::PressureSolver) {
        solver = new PressureSolver();
    } else if (string_solver == SolverTypes::DiffusionTurbSolver) {
        solver = new DiffusionTurbSolver();
    } else if (string_solver == SolverTypes::NSSolver) {
        solver = new NSSolver();
    } else if (string_solver == SolverTypes::NSTurbSolver) {
        solver = new NSTurbSolver();
    } else if (string_solver == SolverTypes::NSTempSolver) {
        solver = new NSTempSolver();
    } else if (string_solver == SolverTypes::NSTempTurbSolver) {
        solver = new NSTempTurbSolver();
    } else if (string_solver == SolverTypes::NSTempConSolver) {
        solver = new NSTempConSolver();
    } else if (string_solver == SolverTypes::NSTempTurbConSolver) {
        solver = new NSTempTurbConSolver();
    } else {
        spdlog::error("Solver not yet implemented! Simulation stopped!");
        std::exit(1);
    }

    // 1. Visualize and test initial conditions
#ifndef PROFILING
    // Solution
    Analysis ana;
    ana.Analyse(solver, 0.);

    // Visualize
    Visual vis;
    vis.Visualize(solver, 0., argv[1]);
#endif

#ifdef _OPENACC
    // Initialize GPU
    acc_device_t dev_type = acc_get_device_type();
    acc_init(dev_type);
#endif

    // 2. Integrate over time and solve numerically
    // Time integration
    TimeIntegration ti(solver, argv[1]);
    ti.run();

    // 3. Compute analytical solution and compare
#ifndef PROFILING
    real t_end = params->getReal("physical_parameters/t_end");
    ana.Analyse(solver, t_end);
#endif

    // 4. Clean up
    delete solver;
    return 0;
}
