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

class FieldIO {
 public:
    explicit FieldIO(FieldController *field_controller);
    void write_out(real t_cur);
    void read(real t_cur, Field *u, Field *v, Field *w, Field *p, Field *T, Field *C);
    void set_filename(std::string &filename) { m_filename = filename; }

 private:
    std::string create_header();
    FieldController *m_field_controller;

    int *m_positions;
    std::string m_filename = "visualisation.dat";
    real m_dt;

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};


#endif /* ARTSS_VISUALISATION_FIELDIO_H */
