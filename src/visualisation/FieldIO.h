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

class FieldIO {
 public:
    FieldIO(FieldController *field_controller, real dt);
    void write_out(real t_cur);
    void read(real t_cur, Field *u, Field *v, Field *w, Field *p, Field *T, Field *C);
    void set_filename(std::string &filename) { m_filename = filename; }

 private:
    std::string create_header();
    FieldController *m_field_controller;

    int m_length;
    std::string m_filename = "visualisation.dat";
    int m_header_length;
    real m_dt;

    int get_position(real t_cur, int length);
};


#endif /* ARTSS_VISUALISATION_FIELDIO_H */
