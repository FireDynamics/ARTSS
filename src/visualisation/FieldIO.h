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
    FieldIO(FieldController &field_controller, real dt);
    void write_out(real t_cur);

 private:
    std::string create_header();

    std::string m_filename = "visualisation.dat";
    int m_header_length;
    FieldController &m_field_controller;
    real m_dt;

    void read(real t_cur);
};


#endif /* ARTSS_VISUALISATION_FIELDIO_H */
