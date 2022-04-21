/// \file       FieldIO.h
/// \brief      Class for reading/writing the raw data of the fields u, v, w, p, T, C
/// \date       Apr 18, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_VISUALISATION_FIELDIO_H
#define ARTSS_VISUALISATION_FIELDIO_H

#include <memory>
#include <string>
#include <iostream>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>

#include "../field/FieldController.h"
#include "../utility/Utility.h"
#include "../solver/SolverController.h"

class FieldIO {
 public:
    explicit FieldIO(const std::string &xml_file_name, const std::string &file_name = "visualisation.dat");
    void write_fields(real t_cur, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C);
    void set_file_name(std::string &file_name) { m_file_name = file_name; }
    void read_fields(real t_cur,
                     const Settings::data_assimilation::field_changes &field_changes,
                     Field &u, Field &v, Field &w,
                     Field &p, Field &T, Field &C);
    void read_fields(const Settings::data_assimilation::field_changes &field_changes,
                     Field &u, Field &v, Field &w,
                     Field &p, Field &T, Field &C);

 private:
    void create_header(HighFive::File &file, const std::string &xml_file_name);
    void read_fields(real t_cur, Field &u, Field &v, Field &w, Field &p, Field &T, Field &C);

    long *m_positions;
    std::string m_file_name;

    long m_pos_time_step;
    std::shared_ptr<spdlog::logger> m_logger;

    void read_field(HighFive::File &file, Field &field);
    void read_vis_field(HighFive::File &file, Field &field, const real t);
};

#endif /* ARTSS_VISUALISATION_FIELDIO_H */
