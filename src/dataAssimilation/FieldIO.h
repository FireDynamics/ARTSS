/// \file       FieldWriter.h
/// \brief      
/// \date       Apr 18, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_VISUALISATION_FIELDIO_H
#define ARTSS_VISUALISATION_FIELDIO_H

#include <string>
#include "../field/FieldController.h"
#include "../utility/Utility.h"
#include "../interfaces/IDataAssimilationWriter.h"
#include "../solver/SolverController.h"

class FieldIO {
 public:
    explicit FieldIO(const SolverController &solver_controller);
    void write(real t_cur, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C);
    void read(std::string &file_name, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C);
    void set_filename(std::string &filename) { m_filename = filename; }
    void read(real t_cur, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C);

 private:
    std::string create_header();

    long *m_positions;
    std::string m_filename = "visualisation.dat";
    real m_dt;
    std::string m_format = ".5e";

    IDataAssimilationWriter *m_func;
    const SolverController &m_solver_controller;
    long m_pos_time_step;
    long m_pos_header;
    int m_length_time_stamp;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

};


#endif /* ARTSS_VISUALISATION_FIELDIO_H */
