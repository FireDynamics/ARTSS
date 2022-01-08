/// \file       FieldIO.h
/// \brief      Class for reading/writing the raw data of the fields u, v, w, p, T, C
/// \date       Apr 18, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_VISUALISATION_FIELDIO_H
#define ARTSS_VISUALISATION_FIELDIO_H

#include <string>
#include "../field/FieldController.h"
#include "../utility/Utility.h"
#include "../solver/SolverController.h"

class FieldIO {
 public:
    FieldIO();
    void write_fields(real t_cur, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C);
    void set_filename(std::string &filename) { m_filename = filename; }
    void read_fields(real t_cur,
                     const Settings::data_assimilation::field_changes &field_changes,
                     Field &u, Field &v, Field &w,
                     Field &p, Field &T, Field &C);

 private:
    std::string create_header();
    void read_fields(real t_cur, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C);

    long *m_positions;
    std::string m_filename = "visualisation.dat";
    std::string m_format = ".5e";

    long m_pos_time_step;
    int m_length_time_stamp;
    std::shared_ptr<spdlog::logger> m_logger;

    void read_field(std::ifstream &file_stream, Field &field);
};


#endif /* ARTSS_VISUALISATION_FIELDIO_H */
